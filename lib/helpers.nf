/*
 * =============================================================================
 * helpers.nf — input parsing and safety checks
 * -----------------------------------------------------------------------------
 * These functions run in the head process, before any task is submitted. They
 * fail early and with a clear message, because a 498-sample run that dies at
 * hour six over a missing manifest wastes a day.
 * =============================================================================
 */

/*
 * Read a boolean parameter that may have come from the command line.
 * -----------------------------------------------------------------------------
 * THIS IS NOT DEFENSIVE PROGRAMMING. It fixes a bug that silently did the
 * opposite of what was asked.
 *
 * A param declared `run_kaiju = true` in nextflow.config is a Boolean. The same
 * param given as `--run_kaiju false` on the command line arrives as the STRING
 * "false" — and every non-empty string is true in Groovy. So `if( params.run_kaiju )`
 * ran Kaiju either way, and the only sign was a task asking for 240 GB of RAM
 * that nobody wanted.
 *
 * Verified on Nextflow 26.04.6: --run_kaiju false yields java.lang.String "false".
 *
 * Use this for EVERY boolean parameter a user can pass, and use it in the config
 * files too, where the same trap applies and helper functions are not available
 * (there, spell it out: params.x.toString().toBoolean()).
 */
def asBool( value ) {
    if( value == null )            return false
    if( value instanceof Boolean ) return value
    return value.toString().trim().toLowerCase() in ['true', 'yes', 'on', '1']
}

/*
 * =============================================================================
 * WHERE THE READS COME FROM — three routes, one answer
 * -----------------------------------------------------------------------------
 * CaMi was written around one ENA study, so its only input route used to be an
 * accession manifest plus a per-run checksum marker. That is the correct route
 * when the reads come from an archive and the wrong route for everybody else:
 * a person who already holds FASTQ — from a sequencing company, or from a
 * download that happened months ago — had to invent accessions to get in.
 *
 * There are now three routes. The mode is decided ONCE, here, and every step
 * downstream sees the same [sample, r1, r2] rows whichever route produced them.
 *
 *   1  --input <file.csv>   a samplesheet. Reads AND, optionally, the design.
 *   2  --input <directory>  a folder of FASTQ, paired by file name.
 *   3  (--input unset)      the ENA manifests. UNCHANGED, and still the default,
 *                           so the study this pipeline was built for keeps
 *                           resolving exactly as it did before.
 *
 * Route 3 is what --step FETCH_BAM, .md5.ok markers and accession-named files
 * belong to. In routes 1 and 2 none of that applies: there is no download step
 * to skip, because there was never a download.
 * =============================================================================
 */

/* 'samplesheet' | 'directory' | 'manifest'. */
def inputMode() {
    if( !params.input )
        return 'manifest'
    def p = file(params.input as String)
    if( !p.exists() )
        error "--input does not exist: ${params.input}"
    return p.isDirectory() ? 'directory' : 'samplesheet'
}

// The optional columns. Absent means "not known", which the statistics step
// reports as untestable rather than guessing at.
def designColumns() {
    return ['individual', 'tissue', 'disease_status', 'dna_prep']
}

/*
 * A sample name becomes part of every output file name, and is compared as a
 * string in every downstream table. Anything a shell would re-interpret, or a
 * name that differs from another only by surrounding space, is rejected here
 * rather than three hours into a run.
 */
def checkSampleName( String name, String where ) {
    if( !name || name != name.trim() )
        error "${where}: sample name is empty or has surrounding space: '${name}'"
    if( !(name ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/) )
        error "${where}: sample name '${name}' is not usable as a file name.\n" +
              "  Use letters, digits, dot, dash and underscore, starting with a letter or digit."
    return name
}

/*
 * Route 1 — the samplesheet.
 * -----------------------------------------------------------------------------
 * Required : sample,fastq_1,fastq_2
 * Optional : individual,tissue,disease_status,dna_prep   (the paired design)
 *
 * Relative FASTQ paths resolve against the SAMPLESHEET's own directory, not the
 * working directory. A sheet that a colleague sends beside their data then works
 * without editing, and the same sheet means the same thing wherever it is run
 * from — which `nextflow run` from an arbitrary directory otherwise breaks.
 *
 * Quoted fields are NOT supported and are rejected with a clear message, and so
 * is any value containing a comma or a tab. A half-correct CSV parser that
 * mis-splits one row in a thousand is worse than none — and a tab-separated
 * sheet may legally hold a comma, which then breaks the design file this writes
 * for the statistics step, a long way from the sheet that caused it.
 *
 * `checkOnly` skips the FASTQ existence check. --step STATS reads nothing from
 * the FASTQ, and requiring them to still be on disk would stop anyone re-running
 * the statistics after archiving the raw reads.
 */
def readSamplesheet( boolean checkOnly = false ) {
    def sheet = file(params.input as String)
    def base  = sheet.parent
    def lines = sheet.readLines().findAll { l -> l?.trim() && !l.trim().startsWith('#') }
    if( lines.size() < 2 )
        error "--input ${sheet} has a header but no rows."

    def sep    = lines[0].contains('\t') ? '\t' : ','
    def header = lines[0].split(sep, -1)*.trim()*.toLowerCase()

    ['sample', 'fastq_1', 'fastq_2'].each { need ->
        if( !header.contains(need) )
            error "--input ${sheet} has no '${need}' column.\n" +
                  "  Required: sample,fastq_1,fastq_2\n" +
                  "  Optional: ${designColumns().join(',')}\n" +
                  "  Found:    ${header.join(',')}"
    }

    def rows = []
    def seen = [] as Set
    lines.drop(1).eachWithIndex { line, i ->
        def lineNo = i + 2
        if( line.contains('"') || line.contains("'") )
            error "--input ${sheet} line ${lineNo}: quoted fields are not supported.\n" +
                  "  Remove the quotes, and keep commas out of paths and labels."
        def f = line.split(sep, -1)*.trim()
        if( f.size() != header.size() )
            error "--input ${sheet} line ${lineNo}: ${f.size()} field(s), header has ${header.size()}."

        def where = "--input ${sheet} line ${lineNo}"

        // A comma in a TAB-separated sheet survives parsing here and then breaks
        // the design file written for --step STATS, with nothing in the error to
        // point back at this line. Reject both separators in every value, so a
        // sheet that parses is a sheet that keeps working downstream.
        def dirty = f.findAll { v -> v.contains(',') || v.contains('\t') }
        if( dirty )
            error "${where}: a value contains a comma or a tab: ${dirty.first()}\n" +
                  "  Neither can be quoted here, and a comma breaks the design file\n" +
                  "  that --step STATS reads. Remove it, or use an underscore."

        def name = checkSampleName( cellAt(header, f, 'sample'), where )
        if( !seen.add(name) )
            error "${where}: sample '${name}' appears more than once.\n" +
                  "  Every sample name must be unique — it names the output files."

        def r1 = resolveRead(header, f, 'fastq_1', base, name, where, checkOnly)
        def r2 = resolveRead(header, f, 'fastq_2', base, name, where, checkOnly)
        if( r1.toString() == r2.toString() )
            error "${where}: fastq_1 and fastq_2 are the same file for '${name}'.\n" +
                  "  CaMi is a paired-end pipeline; single-end input is not supported."

        def design = [:]
        designColumns().each { col -> design[col] = cellAt(header, f, col) }
        rows << [ name, r1, r2, design ]
    }
    return rows
}

/* One cell by column name; '' when the column is absent. */
def cellAt( List header, List fields, String col ) {
    def j = header.indexOf(col)
    return j < 0 ? '' : fields[j]
}

/*
 * One FASTQ path from the sheet, resolved and checked.
 *
 * A path is absolute when it starts with "/" OR names a remote scheme
 * (s3://, gs://, az://, http://, https://, ftp://). Nextflow's file() handles
 * all of those, and resolving "s3://bucket/x.fastq.gz" against the sheet's own
 * directory would produce a local path that cannot exist.
 */
def resolveRead( List header, List fields, String col, base, String name, String where,
                 boolean checkOnly = false ) {
    def raw = cellAt(header, fields, col)
    if( !raw )
        error "${where}: ${col} is empty for sample '${name}'."
    def remote = raw ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/
    def p = ( raw.startsWith('/') || remote ) ? file(raw) : base.resolve(raw)
    // Existence is not checked for a remote URI (it would cost a network call at
    // launch), nor when only the design is wanted — see readSamplesheet.
    if( !checkOnly && !remote && !p.exists() )
        error "${where}: ${col} not found for sample '${name}':\n  ${p}"
    return p
}

/*
 * Route 2 — a directory of FASTQ.
 * -----------------------------------------------------------------------------
 * Recognised mate spellings, which cover what sequencing providers actually
 * send:  <name>_R1_001.fastq.gz  <name>_R1.fastq.gz  <name>_1.fastq.gz
 * and the same four with .fq.gz.
 *
 * A read 1 whose mate is missing is an ERROR, not a skipped file. Half a pair
 * on disk means an interrupted copy or a single-end run, and both deserve to
 * stop the pipeline while somebody is still watching.
 */
def scanFastqDir() {
    def dir = file(params.input as String)
    def mate1 = ~/^(.*?)(_R?1)(_001)?\.(fastq|fq)\.gz$/

    def rows = []
    def bad  = []
    def skipped = []
    dir.list().sort().each { String fname ->
        def m = fname =~ mate1
        if( !m.matches() ) {
            // Anything that LOOKS like a read but was not recognised. Dropping
            // these in silence loses samples while the run still succeeds, which
            // is the failure this whole scan is meant to avoid: uppercase
            // extensions, .fastq without .gz, bzip2, a stray .tmp.
            if( fname ==~ /(?i).*\.f(ast)?q(\..*)?$/ )
                skipped << fname
            return
        }
        def stem = m[0][1]
        def mate = fname.replaceFirst(/(_R?)1(_001)?\.(fastq|fq)\.gz$/) { _all, r, s, ext ->
            "${r}2${s ?: ''}.${ext}.gz"
        }
        def r1 = dir.resolve(fname)
        def r2 = dir.resolve(mate)
        if( !r2.exists() ) { bad << "${fname}  (expected mate ${mate})" ; return }
        rows << [ checkSampleName(stem, "--input ${dir}"), r1, r2, [:] ]
    }

    if( bad )
        error "--input ${dir}: ${bad.size()} read-1 file(s) have no matching read 2:\n" +
              "  " + bad.take(10).join('\n  ') + "\n" +
              "  CaMi is paired-end. Fix the directory, or use a samplesheet with --input <file.csv>."
    if( !rows )
        error "--input ${dir} holds no paired FASTQ.\n" +
              "  Looked for <name>_R1_001, <name>_R1 and <name>_1, with .fastq.gz or .fq.gz." +
              ( skipped ? "\n  ${skipped.size()} file(s) look like reads but do not match: " +
                          skipped.take(5).join(', ') : '' )
    if( skipped )
        log.warn "--input ${dir}: ignoring ${skipped.size()} file(s) that look like reads " +
                 "but do not match a recognised name: ${skipped.take(5).join(', ')}" +
                 ( skipped.size() > 5 ? ", ..." : "" ) +
                 "\n  Recognised: <name>_R1_001, <name>_R1, <name>_1, with .fastq.gz or .fq.gz " +
                 "(lower case). Use a samplesheet to name anything else."

    def dupes = rows*.first().countBy { s -> s }.findAll { _s, n -> n > 1 }.keySet()
    if( dupes )
        error "--input ${dir}: these sample names come from more than one file pair: ${dupes}.\n" +
              "  Two spellings of the same name (e.g. S1_1.fastq.gz and S1_R1.fastq.gz) collide.\n" +
              "  Use a samplesheet with --input <file.csv> to name them apart."
    return rows
}

/* --samples applied to [name, r1, r2, design] rows. */
def selectRows( rows ) {
    if( !params.samples )
        return rows
    def want = wantedSet()
    def kept = rows.findAll { r -> want.contains(r[0]) }
    if( !kept )
        error "--samples matched nothing in --input ${params.input}: ${params.samples}"
    return kept
}

/*
 * Every input row, resolved once: [[sample, r1, r2, accession, source], ...].
 *
 * Reads and provenance come from ONE call, because the samplesheet parse and the
 * directory scan both touch the file system and both log what they found. Two
 * calls meant two scans and a duplicated line in the log.
 *
 * Route 3 rebuilds exactly what it always did — the two disjoint manifests, the
 * collision check between them, and the accession-named paths — because a live
 * 563-sample run resumes against those task hashes.
 */
def resolveInputs() {
    def mode = inputMode()

    if( mode != 'manifest' ) {
        def rows = selectRows( mode == 'samplesheet' ? readSamplesheet() : scanFastqDir() )
        log.info "PROFILE: ${rows.size()} sample(s) from --input (${mode})."
        // No accession exists on these routes. "-" keeps the provenance file one
        // shape whatever the input was, rather than dropping the column.
        return rows.collect { r -> [ r[0], r[1], r[2], '-', mode.toUpperCase() ] }
    }

    def fastqRows = parseSamples()
    def bamRows   = parseBamReadyRows()

    // The two manifests are built from disjoint runs (has-FASTQ vs BAM-only), so
    // their sample codes must never collide. Check anyway: a collision here would
    // silently merge two different animals into one output file.
    def dupes = fastqRows*.first().intersect( bamRows*.first() )
    if( dupes )
        error "Sample code(s) appear in BOTH the FASTQ and BAM-only manifests: ${dupes}.\n" +
              "  This must not happen — check manifests/sample_name_map.csv and ${params.bam_manifest}."

    log.info "PROFILE: ${fastqRows.size()} FASTQ sample(s) + ${bamRows.size()} BAM-derived sample(s)."

    def rows = fastqRows.collect { sample, acc ->
        def r1 = file("${params.fastq_dir}/${acc}_1.fastq.gz")
        def r2 = file("${params.fastq_dir}/${acc}_2.fastq.gz")
        if( !r1.exists() || !r2.exists() )
            error "Missing FASTQ for ${sample} (${acc}) in ${params.fastq_dir}"
        [ sample, r1, r2, acc, 'FASTQ' ]
    }
    rows += bamRows.collect { sample, acc ->
        [ sample,
          file("${params.bamdir}/${sample}_1.fastq.gz"),
          file("${params.bamdir}/${sample}_2.fastq.gz"),
          acc, 'BAM' ]
    }
    return rows
}

/*
 * The design columns, for --step STATS: [[sample, individual, tissue,
 * disease_status, dna_prep], ...], or null when this run has no samplesheet.
 *
 * Null is not "no design". It means route 3, where the design comes from the ENA
 * manifests and the sample-code grammar exactly as it always has. Returning an
 * empty table instead would quietly replace a working design with nothing.
 */
def resolveDesign() {
    if( inputMode() == 'manifest' )
        return null
    // checkOnly: --step STATS reads the classifier output, never the FASTQ, so
    // the reads may have been archived since PROFILE ran.
    def rows = selectRows( inputMode() == 'samplesheet' ? readSamplesheet(true) : scanFastqDir() )
    return rows.collect { r -> [ r[0] ] + designColumns().collect { c -> (r[3][c] ?: '') } }
}

/*
 * Read manifests/sample_name_map.csv -> [[sample_code, run_accession], ...].
 * The CSV has one row per mate, so it collapses to one row per run.
 * Applies --samples (codes or accessions) and the .md5.ok requirement.
 */
def parseSamples() {
    def mapFile = file(params.sample_map)
    if( !mapFile.exists() )
        error "Sample map not found: ${params.sample_map}"

    def byAcc = [:]
    mapFile.readLines().drop(1).each { line ->
        if( !line?.trim() ) return
        def f = line.split(',')
        if( f.size() < 2 ) return
        byAcc[ f[0].trim() ] = f[1].trim()
    }
    if( !byAcc )
        error "No rows parsed from ${params.sample_map}"

    if( params.samples ) {
        def want = wantedSet()
        byAcc = byAcc.findAll { acc, code -> want.contains(code) || want.contains(acc) }
        if( !byAcc )
            error "--samples matched nothing in ${params.sample_map}: ${params.samples}"
    }

    def rows = []
    def notReady = 0
    byAcc.each { acc, code ->
        // A file of the right size can still be corrupt (README trap 2), so the
        // MD5 marker is required, not the file itself.
        if( asBool(params.require_md5ok) && !file("${params.fastq_dir}/${acc}.md5.ok").exists() ) {
            notReady += 1
            return
        }
        rows << [ code, acc ]
    }

    if( notReady > 0 )
        log.info "Skipped ${notReady} run(s): no verified .md5.ok marker yet."
    if( !rows )
        error "No runs with a verified .md5.ok marker in ${params.fastq_dir}.\n" +
              "  The download may still be in progress. Check with scripts/download/03_status.sh."
    return rows.sort { r -> r[0] }
}

/* Read manifests/skipped_bam_only.tsv -> [[sample_code, run_accession], ...]. */
def parseBamRuns() {
    def f = file(params.bam_manifest)
    // A missing manifest means "no archive run needs the BAM route", which is the
    // ordinary case for every cohort that is not PRJEB58149. It used to be a hard
    // error, so `--step all` died immediately on anyone else's data. Absence is
    // reported and returns nothing, the same way parseBamReadyRows() treats it.
    if( !f.exists() ) {
        log.info "FETCH_BAM: no BAM manifest at ${params.bam_manifest} — nothing to fetch."
        return []
    }

    def rows = []
    f.readLines().drop(1).each { line ->
        if( !line?.trim() ) return
        def c = line.split('\t')
        if( c.size() < 2 ) return
        rows << [ c[1].trim(), c[0].trim() ]
    }

    // --samples is lenient here on purpose. The BAM-only runs and the FASTQ runs
    // are disjoint sets, so a selection aimed at PROFILE normally matches no BAM
    // run at all. That is not an error — it just means there is nothing to fetch.
    // --step FETCH_BAM checks for an empty selection itself.
    if( params.samples ) {
        def want = wantedSet()
        rows = rows.findAll { code, acc -> want.contains(code) || want.contains(acc) }
    }

    // Idempotent: skip runs already fetched.
    def todo = rows.findAll { code, _acc -> !file("${params.bamdir}/${code}_1.fastq.gz").exists() }
    log.info "FETCH_BAM: ${todo.size()} run(s) to fetch, ${rows.size() - todo.size()} already present."
    return todo
}

/*
 * The other half of parseBamRuns(): BAM-only samples that HAVE already been
 * fetched, so PROFILE can pick them up. Without this, FETCH_BAM writes FASTQ
 * into fastq_from_bam/ and nothing downstream ever reads it — the 65 runs
 * (64 of them tumour) would be fetched and then silently ignored.
 * --samples applies the same way it does in parseSamples().
 */
def parseBamReadyRows() {
    def f = file(params.bam_manifest)
    if( !f.exists() )
        return []   // no manifest yet: nothing to add, not an error

    def rows = []
    f.readLines().drop(1).each { line ->
        if( !line?.trim() ) return
        def c = line.split('\t')
        if( c.size() < 2 ) return
        rows << [ c[1].trim(), c[0].trim() ]
    }

    if( params.samples ) {
        def want = wantedSet()
        rows = rows.findAll { code, acc -> want.contains(code) || want.contains(acc) }
    }

    def ready = rows.findAll { code, _acc ->
        file("${params.bamdir}/${code}_1.fastq.gz").exists() && file("${params.bamdir}/${code}_2.fastq.gz").exists()
    }
    if( ready )
        log.info "PROFILE: ${ready.size()} additional BAM-derived sample(s) from ${params.bamdir}."
    return ready
}

/* --samples is either a comma-separated list or a file with one entry per line. */
def wantedSet() {
    def s = params.samples.toString()
    return ( file(s).exists() ? file(s).readLines()*.trim().findAll{ v -> v }
                              : s.split(',')*.trim() ) as Set
}

/*
 * Fail before submitting anything, rather than at hour six of a 563-sample run.
 * -----------------------------------------------------------------------------
 * Everything here runs in the head process, so the cost is milliseconds and the
 * error arrives while the person who caused it is still watching.
 */
def checkLocation() {
    checkScratch()
    checkInputs()
    checkKaijuInputs()
}

/*
 * There must be a way in, and the error must name every way in.
 *
 * The old message pointed only at manifests/sample_name_map.csv, which reads as
 * "this pipeline needs an ENA study" to somebody who simply has FASTQ on disk.
 */
def checkInputs() {
    def scheme = params.sample_code_scheme
    if( scheme && !(scheme in ['bruzos', 'none']) )
        error "--sample_code_scheme must be 'bruzos' or 'none', not '${scheme}'."

    if( !(params.step in ['PROFILE', 'all']) )
        return
    if( params.input ) {
        inputMode()          // fails here if the path does not exist
        return
    }
    if( file(params.sample_map).exists() )
        return

    error """
    No reads to work on. CaMi takes them by one of three routes:

      1. A samplesheet — reads, and optionally the paired design:
             --input samples.csv
         with the columns  sample,fastq_1,fastq_2
         and, when you want --step STATS,  ${designColumns().join(',')}

      2. A directory of FASTQ, paired by file name:
             --input /path/to/fastq
         Recognised: <name>_R1_001, <name>_R1, <name>_1, .fastq.gz or .fq.gz

      3. An ENA study, downloaded first with scripts/download/:
         ${params.sample_map}
         which does not exist.

    Route 1 or 2 needs no download and no manifest. See docs/usage.md.
    """.stripIndent()
}

/*
 * OFF unless a site profile sets --enforce_scratch true.
 *
 * It exists for a cluster whose data sits on a node's LOCAL disk, where the same
 * bytes have two names: /scratch/<user>/... on the node itself, and something
 * like /mnt/<node>/<user>/... from the login node. A head process started on the
 * wrong side hands every task a path that does not exist there, so every task
 * fails; and because -resume keys on the path, the two spellings would also split
 * the cache in half and repeat days of finished work.
 *
 * On ordinary shared storage none of that applies, which is why the default is
 * false. See conf/site/mpi_bremen.config for the site that needs it.
 */
def checkScratch() {
    if( !asBool(params.enforce_scratch) )
        return
    def root = file(params.projectRoot).toAbsolutePath().normalize().toString()
    if( !root.startsWith('/scratch/') )
        error """
        --enforce_scratch is set, but the data root is ${root}, not under /scratch.

        This site keeps the data on a compute node's local disk. The Nextflow head
        process must run ON that node, so that it and its tasks spell the path the
        same way. Submit it with:
            ./run_pipeline.sh --step <NAME>
        """.stripIndent()
}

/*
 * Kaiju needs three paths and gives an unhelpful error for any missing one, deep
 * inside a task, after the whole host-removal chain has already run for hours.
 * The index is ~187 GB and is never bundled, so a fresh installation reaching
 * --step PROFILE with kaiju_db unset is the expected mistake, not an exotic one.
 */
def checkKaijuInputs() {
    if( !asBool(params.run_kaiju) || !(params.step in ['PROFILE', 'all']) )
        return
    def missing = ['kaiju_db', 'kaiju_nodes', 'kaiju_names'].findAll { k ->
        !params[k] || !file(params[k] as String).exists()
    }
    if( missing )
        error """
        --run_kaiju is true, but these are unset or do not exist: ${missing.join(', ')}

        Kaiju's nr_euk index is about 187 GB, so CaMi never downloads or bundles
        it. Point the pipeline at your site's copy:
            --kaiju_db    /db/Kaiju/kaiju_db_nr_euk.fmi
            --kaiju_nodes /db/Kaiju/nodes.dmp
            --kaiju_names /db/Kaiju/names.dmp
        or turn the step off with --run_kaiju false. docs/usage.md says how to
        build the index if your site has no copy.

        Turning it off is not free: on identical reads Kaiju found 3,249 bacterial
        pairs where Kraken 2 found 210. See modules/local/kaiju.nf.
        """.stripIndent()
}

/* The bwa-mem2 host index, as built by the REFERENCE step. */
def existingHostIdx() {
    def fs = files("${params.refdir}/${params.host_prefix}.{0123,amb,ann,bwt.2bit.64,pac}")
    if( fs.size() < 5 )
        error "Host index incomplete in ${params.refdir} (found ${fs.size()} of 5 files).\n" +
              "  Build it first:  nextflow run ${workflow.projectDir} --step REFERENCE\n" +
              "  (or ./run_pipeline.sh --step REFERENCE if your site has a wrapper)"
    return fs
}
