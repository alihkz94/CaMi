/*
 * BUILD_MM2_INDEX — minimap2 index, used for the human genome and the vibrio panel
 */
process BUILD_MM2_INDEX {
    tag "$name"

    input:
    tuple val(name), path(fagz)

    output:
    path "${name}.mmi"

    script:
    """
    minimap2 -x sr -t ${task.cpus} -d ${name}.mmi ${fagz}
    """
}
