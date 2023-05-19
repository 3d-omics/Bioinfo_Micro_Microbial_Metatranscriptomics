def get_star_out_prefix(wildcards):
    return STAR / f"{wildcards.sample}.{wildcards.library}."


def get_star_output_r1(wildcards):
    return STAR / f"{wildcards.sample}.{wildcards.library}.Unmapped.out.mate1"


def get_star_output_r2(wildcards):
    return STAR / f"{wildcards.sample}.{wildcards.library}.Unmapped.out.mate2"
