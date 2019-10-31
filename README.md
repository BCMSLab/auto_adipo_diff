# auto_adipo_diff

Code for the Cells paper (Transcriptional Regulation of Autophagy Genes via Stage-Specific Activation of CEBPB and PPARG
during Adipogenesis)

## Setting up the docker environment

The analysis was run on a [docker](https://hub.docker.com/r/bcmslab/autoreg/) image based on the the latest **rocker/verse**.
Other R packages were added to the image and were made available as an image that can be obtained and launched on any local
machine running [docker](https://hub.docker.com/r/bcmslab/autoreg/).

```bash
$ docker pull bcmslab/autoreg:latest
$ docker run -it bcmslab/autoreg:latest bash
```

## Obtaining the source code

The source code is hosted publicly on this repository in the form of a research compendium. This includes the scripts to 
reproduce the figures and tables in this manuscript. Another public repository [autoreg](https://github.com/BCMSLab/autoreg)
contains the scripts to download, prepare and analyze the data. From within the container, [git](https://git-scm.com) can be
used to clone the source code using the `--recurse-submodules` flag to include the submodules.

The following code clones the repository containing the source code.

```bash
$ git clone --recurse-submodules http://github.com/BCMSLab/auto_adipo_diff
```

## Runing the analysis

In the submodule `autoreg`, run `make`

```bash
$ cd auto_adipo_diff/autoreg
$ make analysis
```

## Generating figures and tables

The script to generate the figures and tables in the manuscirpt can be als ran through `make`

```bash
$ cd auto_adipo_diff
$ make all
```

## Details of the R environment
The version of **R** that was used to perform this analysis is the 3.7.0 (2019-04-07) on `x86\_64-pc-linux-gnu`.

## More

This manuscript was published under the title [Transcriptional Regulation of Autophagy Genes via Stage-Specific 
Activation of CEBPB and PPARG during Adipogenesis: A Systematic Study Using Public Gene Expression and Transcription 
Factor Binding Datasets](https://www.mdpi.com/2073-4409/8/11/1321)
