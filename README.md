# Overview

Ptolemy is a reference-free approach for analysing microbial genome architectures. In a nutshell, it uses a "top-down" approach to align multiple genomes via synteny analysis of gene annotations. It requires a set of FASTA-formatted-assemblies and corresponding GFF-formatted-annotations. The resulting alignment is represented as a graph describing structural similarities and differences of (sub-)populations.

Ptolemy has been accepted as a conference proceeding paper at [ECCB 2018](http://eccb18.org/proceedings/) and will be published in [Bioinformatics](https://academic.oup.com/bioinformatics).

# Executable JAR

Executable jar files are available under [releases](https://github.com/AbeelLab/ptolemy/releases/latest).

# Running Ptolemy

Ptolemy requires a tab-delimited file containing unique sample identifier, path to misc.assembly, and path to gene annotations. For example:

```
Genome1 path/to/misc.assembly/genome1.fa path/to/annotations/genome1.gff
Genome2 path/to/misc.assembly/genome2.fa path/to/annotations/genome2.gff
Genome3 path/to/misc.assembly/genome3.fa path/to/annotations/genome3.gff
```

There are then three main steps in Ptolemy:
1. Database creation ( *java -jar ptolemy.jar extract ...* )
2. Multiple-genome alignment via syntenic anchoring ( *java -jar ptolemy.jar syntenic-anchors ...* )
3. Canonical graph construction ( *java -jar ptolemy.jar canonical-quiver ...* )

A typical workflow:

```
java -jar ptolemy.jar extract -g genome_list.txt -o ptolemy_db
java -jar ptolemy.jar syntenic-anchors --db ptolemy_db -o  .
java -jar ptolemy.jar canonical-quiver -s syntenic_anchors.txt --db ptolemy_db -o .
```

The graph is stored as a [GFA-formatted file](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) and can be visualized via graph-visualizers such as [Bandage](https://rrwick.github.io/Bandage/).

Test-data available under 'testing_data' directory which contains full [Pacbio assemblies](https://yjx1217.github.io/Yeast_PacBio_2016/data/) of a single yeast chromosome from three genomes.
