/*
 * FETCH_GENOME — download a reference genome
 * -----------------------------------------------------------------------------
 * aria2c uses many connections. The site proxy flaps, so the retry settings in
 * conf/base.config matter more than the speed here.
 */
process FETCH_GENOME {
    tag "$name"

    input:
    tuple val(name), val(url)

    output:
    tuple val(name), path(name)

    script:
    """
    aria2c -x16 -s16 -c --console-log-level=warn --auto-file-renaming=false \\
           -o ${name} "${url}"
    """
}
