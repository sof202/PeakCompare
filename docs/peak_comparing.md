# Peak Comparing

After completing the [peak calling pipeline](./peak_calling.md) for both your
reference dataset and your comparison dataset, you can now compare regions
between the two. The purpose of this script is to generate a metric that aims
to tell the user how likely that a comparable peak in the comparison dataset
exists (given a peak exists in the reference dataset). Details on how to
interpret the output metric can be found [here](./metric_interpretation.md).

## Running

Before running this script, it is assumed you have ran the peak calling 
pipeline for both of your datasets. After accomplishing this, identify a peak
(or selection of close by peaks) in the reference dataset (noting down the 
chromosome, start and end base positions). Feed this information into the
configuration file `example_config_compare.txt`. This configuration file asks
for a few different parameters, head over to [this section] to understand what
these parameters mean and how they are used.

After filling out the configuration file, submit the job to a SLURM queue with

```bash
sbatch -p queue .../PeakCompare.sh path/to/config_file.txt
```

### Ignoring this script

This script is only given to provide a means for the user to submit the job
into a SLURM queue. You could alternatively just run the underlying python
script yourself. To do this, activate the conda environment with:

```bash
conda activate PeakCompare-MACS
```

Then run the python script called `peak_compare.py` which can be found in the
`Python_Scripts` directory of this repository. There are quite a few arguments
that this script requires due to the multitude of files it needs to read in
(another reason behind the wrapper script). For this reason, you may decide to
create your own wrapper script that achieves your own goals (perhaps you want
to look at multiple regions in the genome simultaneously).

#### grep

To save time on reading in files, the wrapper script uses grep to pull out the
chromosome of interest before passing them to `peak_compare.py`. If you plan on
running the python script interactively (or you are writing your own wrapper
script), consider implementing this as well.

## Parallelism

The python script here allows for inspecting multiple regions at once. Give
a comma separated list of chromosomes, start bases and end bases for each 
region you want to inspect. The order of the results may not be in the same
order as your inputs. Because of this, the output of the script is 4 tab
separated columns.

|Chromosome|Start|End|Metric|
|----------|-----|---|------|

This is written directly to the standard output stream, so you can use basic
linux commands to work with it (awk, grep, `>`, *etc.*).

## How the metric is calculated

The metric is a simple ratio of the number of bases in psuedopeaks in the
comparison and the number of bases in peaks in the reference dataset.
Pseudopeaks is a term that describes regions in your comparison set that are
'almost' peaks when compared to the reference dataset. A more formal definition
is given [here](#pseudopeak-definition).

As alluded to in [cutoff analysis](./peak_calling.md#cutoff-analysis), cutoffs
are tricky little things. No matter what cutoff you choose, you will almost
always be classing some regions as peaks because they only just exceed the
cutoff (and classing some regions not as peaks for the opposite reason). This
is problematic when comparing two datasets. How can we determine if a peak is
truly absent in the comparison dataset when there is a possibility that a
single extra read in the region could tip the balance? 

The solution given with this tool is to look at the confidence intervals for
each p-value instead of looking at whether they hit certain thresholds. If the
confidence intervals for each dataset overlap, it suggests that (under the
chosen [significance](#significance)) the true value could be shared by each
dataset at this point. Considering these intervals is both more lax and strict
than simply using a single threshold value. If the region in the reference
dataset is a very significant peak (very low p-value), then this method
required the same region in the comparison dataset to have similarly small
p-values to be considered as a pseudopeak.

Comparing confidence intervals is useful, but comes with a drawback in this
particular case of using MACS to determine peaks. As mentioned in 
[this section](./peak_calling.md#why-merged-and-unmerged-peaks), some peaks
may have a very high p-value (not significant), and yet be called a peak due to
the merging process. This is why both files are required by this python script
and allows us to now formally state what defines a psuedopeak.

### Pseudopeak definition

Def<sup>n</sup>: A base position is defined to be in a pseudopeak region if
it satisfies either of these criteria:

(i) The base position is within a peak region in the unmerged reference dataset
peaks, and the confidence interval in the comparison dataset overlaps or
exceeds the reference dataset's confidence interval (for this base position).

(ii) The base position surpasses the cutoff threshold in the comparison dataset
and is only within a peak region in the merged reference dataset peaks (not
present in the unmerged peaks).

### Confidence interval calculation

It should be noted first that p-values from MACS are calculated using the
Poisson distribution (as pileup counts generally follow this distribution).
Now, the confidence interval calculation is carried out in two steps:

1) A confidence interval is calculated for each position in the bias track.
This requires the following steps:
  - The local variance (see [window size](#window-size)) for the bias track is
  calculated.
  - The local standard error is calculated.
  - The confidence interval is found by using $\hat\lambda = \lambda \pm
  z_\alpha \cdot std\_err$. Where $z_\alpha$ is the critical value for the
  given [significance](#significance).
2) A confidence interval for the p-value is calculated for each position using
the Poisson distribution with arguments:
  - The $\hat\lambda$ confidence interval calculated in 1)
  - The number of reads in that position seen in the coverage (pileup) track

## Parameters

The underlying python script requires a few parameters that can be specified
in the configuration file. These are:

- [The cutoff](#cutoff)
- [The significance](#significance)
- [The window size](#window-size)
- ['unmerged'](#unmerged)

### Cutoff

This is the cutoff that was used with the reference dataset when peak calling.
If you didn't set this value yourself (or simply forgot), you can find out
what cutoff value was used by looking at the log file for the job you
submitted. This should look like:

```
───────┬───────────────────────────────────
       │ File: 01-Jan~00-12345678.log
───────┼───────────────────────────────────
   1   │ Using a cutoff value of 5.
───────┴───────────────────────────────────
```

### Significance

This is the significance used when calculating the confidence intervals
described [here](#how-the-metric-is-calculated). Picking a value that is
smaller will result in tighter confidence intervals and therefore results in
stricter [criteria for psuedopeaks](#pseudopeak-definition). There is no
recommended value here, just be aware that confidence intervals are being
compared against each other (so a significance of 95% is more like 90%).

### Window size

In calculating the confidence intervals used to determine pseudopeaks, a window
is used to calculate local variance in the bias track. A suitable value for
the window size will be between a quarter and half of the fragment length for
each dataset. It should be noted that a smaller window size will result in the
variance calculated being more localised. However, this comes with the caveat
that less variance can be observed (consider the case when window size is 0). A
careful balance should be chosen here.

### Unmerged

You may wish to discount peak regions in the reference dataset that are a
result of the merging process entirely when calculating the metric. To do this
you can specify `UNMERGED=1` in the configuration file (or use `--unmerged` if
using the python script). This will mean that the metric will only look at the
ratio of psuedopeaks (in the comparison dataset) with the unmerged peaks (in
the reference dataset). All peaks that are a result of merging will be ignored.

## Time

Running this script is fast (a couple of seconds). The main slowdown comes with
reading in files. This is why [grep](#grep) is utilised in the provided wrapper
script. The complexity of your data can also increase read times as bedgraph
files (output of most MACS commands) can massively vary in size due to this
factor. The bedgraph file at minimum can be 23 lines long (one for each
chromosome) up to ~2,700,000,000 lines long (one for each base pair).
