**Prerequisites**
=================

## Docker

The computational pipeline is contained in a docker image. We strongly recommend using it for your analyses. Follow the official docs for installation instructions: [MacOS](https://docs.docker.com/docker-for-mac/install/) or [Linux](https://docs.docker.com/install/linux/docker-ce/ubuntu/). Unfortunately, our provided pipeline does not support Windows.

Next, clone this repository and build the docker image.

```
git clone https://github.com/octantbio/octopus3.git
cd octopus3
docker build -t octopus3:release .
```

It is possible to run OCTOPUS without `docker`, but not recommended. If you insist, see the [Dockerfile](/Dockerfile) for help.

## Hardware

For optimal performance, we recommend deploying on a machine with at least 32 GB of RAM and 16 cores. We use [GNU Parallel](https://www.gnu.org/software/parallel/) to distribute tasks where possible, and empirically, the RAM requirements seems to scale with the number of cores (e.g. 64 GB of RAM for 64 cores). The compute times (for 384 wells) scale asymptotically, suggesting a potential disk bottleneck.

```
64 cores -> 64 GB RAM -> 13 Mins
32 cores -> 32 GB RAM -> 15 Mins
24 cores -> 24 GB RAM -> 20 Mins
16 cores -> 16 GB RAM -> 25 Mins
8  cores ->  8 GB RAM -> 45 Mins
```

We performed all trials on two Intel(R) Xeon(R) Gold 6130 CPUs @ 2.10GHz, limiting the cores through the docker `--cpuset-cpus=N` flag, and estimated peak RAM usage with `/bin/free`.

