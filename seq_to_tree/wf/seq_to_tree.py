"""
align sequences using MAFFT
"""

from enum import Enum
import os
from pathlib import Path
import subprocess
from typing import Optional

from latch import small_task, workflow
from latch.types import LatchFile, LatchDir


class AlignmentMode(Enum):
    linsi = "L-INS-i"
    fftns2 = "FFT-NS-2"
    auto = "auto"

class TrimmingMode(Enum):
    gappy = "gappy"
    smart_gap = "smart-gap"
    kpi = "kpi"
    kpi_gappy = "kpi-gappy"
    kpi_smart_gap = "kpi-smart-gap"
    kpic = "kpic"
    kpic_gappy = "kpic-gappy"
    kpic_smart_gap = "kpic-smart-gap"


@small_task
def align_sequences_task(
    unaligned_seqs: LatchFile,
    alignment_mode: AlignmentMode = AlignmentMode.auto,
    output_prefix: Optional[str] = "seq2tree"
    ) -> LatchFile:

    ## specify output file name
    alignment_out_file = Path(f"{output_prefix}_alignment_mafft.fa").resolve()

    #local directory to put output files in
    local_dir = "/root/seq_to_tree/" 

    ## logic for how to align seqs
    if alignment_mode.value == AlignmentMode.linsi:
        _mafft_cmd = [
            "mafft-linsi",
            unaligned_seqs.local_path,
        ]
    elif alignment_mode.value == AlignmentMode.fftns2:
        _mafft_cmd = [
            "mafft",
            unaligned_seqs.local_path,
        ]
    else:  # alignment_mode.value == AlignmentMode.auto
        _mafft_cmd = [
            "mafft",
            "--auto",
            unaligned_seqs.local_path,
        ]

    with open(alignment_out_file, "w") as f:
        subprocess.call(_mafft_cmd, stdout=f)

    return LatchFile(str(alignment_out_file), f"latch:///{local_dir}{alignment_out_file.name}")


@small_task
def trim_alignment_task(
    alignment: LatchFile,
    output_prefix: Optional[str] = "seq2tree",
    gap_threshold: float = 0.9,
    trimming_mode: TrimmingMode = TrimmingMode.smart_gap
    ) -> LatchFile:

    # Set output file name
    trimmed_aln_fasta = Path(f"{output_prefix}_trimmed_clipkit.fa").resolve()

    #local directory to put output files in
    local_dir = "/root/seq_to_tree/" 

    # Set gap threshold is user did not provide a value
    if not gap_threshold:
        gap_threshold=0.9

    _clipkit_cmd = [
        "clipkit",
        alignment.local_path,
        "-o",
        str(trimmed_aln_fasta),
        "-m",
        trimming_mode.value,
        "-g",
        str(gap_threshold)
    ]

    subprocess.run(_clipkit_cmd)

    return LatchFile(str(trimmed_aln_fasta), f"latch:///{local_dir}{trimmed_aln_fasta.name}")


@small_task
def infer_phylogeny_task(
    trimmed_alignment: LatchFile,
    output_dir: LatchDir,
    output_prefix: Optional[str] = "seq2tree",
    ufboot_reps: Optional[int] = 1000
    ) -> LatchDir:

    # if no output prefix is specified, assign output prefix to "seq2tree"
    if isinstance(output_prefix, type(None)):
        output_prefix = "seq2tree"

    #local directory to put output files in
    local_dir = "/root/seq_to_tree/" 
    # iqtree prefix including local path
    local_prefix = os.path.join(local_dir, output_prefix) 

    # iqtree command
    _iqtree_cmd = [
        "iqtree2",
        "-s",
        trimmed_alignment.local_path,
        "-pre",
        str(local_prefix),
        "-nt",
        "AUTO",
        "-m",
        "TEST",
        "-bb",
        str(ufboot_reps)
    ]

    subprocess.run(_iqtree_cmd)
    return LatchDir(local_dir, output_dir.remote_path) 


@workflow
def seq_to_tree(
    unaligned_seqs: LatchFile,
    output_dir: LatchDir,
    alignment_mode: AlignmentMode = AlignmentMode.auto,
    trimming_mode: TrimmingMode = TrimmingMode.smart_gap, 
    gap_threshold: float = 0.9,
    ufboot_reps: Optional[int] = 1000,
    output_prefix: Optional[str] = "seq_to_tree"
    ) -> LatchDir:
    """
    Seq to tree
    ----
    # Seq_to_tree
    ## About
    The Seq_to_tree workflow goes from a multi-FASTA file of sequences to a phylogenetic tree. 
    Specifically, the seq_to_tree workflow will align sequences using Mafft, trim the alignment
    using ClipKIT, and then infer the evolutionary history of the sequences using IQTREE.

    Users can modify many, but not all, components of running each software. These are described
    in detail below. However, user's must name their output directory and input their multi-FASTA
    file. Otherwise, default parameters will be selected for all other options.

    <br /><br />

    If you found seq_to_tree useful, please cite *MAFFT Multiple Sequence Alignment
    Software Version 7: Improvements in Performance and Usability*. Katoh & Standley 2013,
    Molecular Biology and Evolution. doi:
    [10.1093/molbev/mst010](https://academic.oup.com/mbe/article/30/4/772/1073398)*;
    *ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference.
    Steenwyk et al. 2020, PLoS Biology. doi: 
    [10.1371/journal.pbio.3001007](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007)*;
    and *IQ-TREE 2: New models and efficient methods for phylogenetic inference in the
    genomic era. Minh et al. 2020, Molecular Biology and Evolution. doi:
    [10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015)*.

    <br /><br />

    ## Mafft options
    Herein, the various aligning modes implemented in MAFFT are described. If you are not sure
    which is appropriate for you, I recommend using the auto trimming mode.
    <br />
    - L-INS-i: an accurate option (L-INS-i) for an alignment of up to ∼200 sequences × ∼2,000 sites
    - FFT-NS-2: a fast option (FFT-NS-2) for a larger sequence alignment
    - auto: allow MAFFT to decide

    ## ClipKIT options
    ClipKIT can be run with eight different modes. If you are not sure which is appropriate for you,
    I recommend trimming using the smart-gap mode.
    <br />
    smart-gap: dynamic determination of gaps threshold
    gappy: trim all sites that are above a threshold of gappyness (default: 0.9)
    kpic: keep only parismony informative and constant sites
    kpic-smart-gap: a combination of kpic- and smart-gap-based trimming
    kpic-gappy: a combination of kpic- and gappy-based trimming
    kpi: keep only parsimony informative sites
    kpi-smart-gap: a combination of kpi- and smart-gap-based trimming
    kpi-gappy: a combination of kpi- and gappy-based trimming

    ## IQTREE2 options
    The majority of parameters for IQTREE2 are automatically determined (e.g., substitution model).
    The user can flexibly define how many UFBoot replicates to use when assessing bipartition support.
    <br />

    __metadata__:
        display_name: Align sequences using MAFFT
        author: Jacob L. Steenwyk
            name: Jacob L. Steenwyk
            email: jlsteenwyk@gmail.com
            github: https://github.com/JLSteenwyk
        repository: https://mafft.cbrc.jp/alignment/software/
        license:
            id: BSD

    Args:

        output_dir:
            Output directory
			__metadata__:
				display_name: "Output directory"
				appearance:
					comment: "Output directory"
        
        output_prefix:
            Prefix of all outputted files.
			__metadata__:
				display_name: "Prefix of outputted files."
				appearance:
					comment: "Prefix of outputted files."

        unaligned_seqs:
            Input multi-FASTA file of nucleotide or amino acid sequences 
            __metadata__:
                display_name: "Input multi-FASTA file"
                appearance:
					comment: "Input multi-FASTA file"

        alignment_mode:
            Mode for multiple sequence alignment 
            __metadata__:
                display_name: "Alignment mode (see About for details)"
                appearance:
					comment: "- auto
                    - L-INS-i
                    - FFT-NS-2"

        trimming_mode:
            Mode used for trimming. See "About" for more information.
            __metadata__:
                display_name: "Trimming mode"
                appearance:
                    comment: "- smart-gap: dynamic determination of gaps threshold
                        - gappy: trim all sites that are above a threshold of gappyness (default: 0.9)
                        - kpic: keep only parismony informative and constant sites
                        - kpic-smart-gap: a combination of kpic- and smart-gap-based trimming
                        - kpic-gappy: a combination of kpic- and gappy-based trimming
                        - kpi: keep only parsimony informative sites
                        - kpi-smart-gap: a combination of kpi- and smart-gap-based trimming
                        - kpi-gappy: a combination of kpi- and gappy-based trimming"

        gap_threshold:
            Specifies gaps threshold (default: 0.9). Ignored if smart-gap is used.
			__metadata__:
				display_name: "Gappyness threshold"
				appearance:
					comment: "Specifies gaps threshold (default: 0.9). Ignored if smart-gap is used."

        ufboot_reps:
            UFBoot (Ultrafast bootstrap approximation replicates)
            __metadata__:
                display_name: "Integer for number of UFBoot replicates to run"
                appearance:
					comment: "UFBoot replicates"

    """

    alignment = align_sequences_task(
        unaligned_seqs=unaligned_seqs,
        alignment_mode=alignment_mode,
        output_prefix=output_prefix
        )

    trimmed_alignment = trim_alignment_task(
        alignment=alignment,
        output_prefix=output_prefix,
        gap_threshold=gap_threshold,
        trimming_mode=trimming_mode
    )

    return infer_phylogeny_task(
        trimmed_alignment=trimmed_alignment,
        output_dir=output_dir,
        output_prefix=output_prefix,
        ufboot_reps=ufboot_reps
    )
    