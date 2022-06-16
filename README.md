<p align="center">
    <a href="https://latch.bio/">
        <img src="images/latch_logo.png" width=37.5%>
    </a>
    <br />
    <font width=15%>Mafft -> ClipKIT -> IQTree2</font>
    <br /><br />
    <a href="https://console.latch.bio/explore/62302/info">
    <span style="font-size:larger;">Click here to see the workflow!</span></a>
    </br></br>
    Workflow author: <a href="https://jlsteenwyk.com/">Jacob L. Steenwyk</a>
    </br>
    <a href="https://twitter.com/intent/follow?screen_name=jlsteenwyk" alt="Author Twitter">
        <img src="https://img.shields.io/twitter/follow/jlsteenwyk?style=social&logo=twitter"
            alt="follow on Twitter">
    </a>
</p>

</br>

>This is an implementation of a 'typical' workflow for inferring the evolutionary history among a set of sequences.

</br>

# Seq_to_tree
## About
The Seq_to_tree workflow goes from a multi-FASTA file of sequences to a phylogenetic tree. Specifically, the seq_to_tree workflow will align sequences using Mafft, trim the alignment using ClipKIT, and then infer the evolutionary history of the sequences using IQTREE.

Users can modify many, but not all, components of running each software. These are described in detail below. However, user’s must name their output directory and input their multi-FASTA file. Otherwise, default parameters will be selected for all other options.

## Citations
If you found seq_to_tree useful, please cite 
- MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. Katoh & Standley 2013, Molecular Biology and Evolution. doi: [10.1093/molbev/mst010](https://academic.oup.com/mbe/article/30/4/772/1073398)*
- ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. Steenwyk et al. 2020, PLoS Biology. doi: [10.1371/journal.pbio.3001007](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007)
- IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Minh et al. 2020, Molecular Biology and Evolution. doi: [10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015).