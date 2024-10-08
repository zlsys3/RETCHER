% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.R
\name{runSampleGroup}
\alias{runSampleGroup}
\title{Run Samples Grouping}
\usage{
runSampleGroup(
  ccfInputList,
  sampleName = NULL,
  minPresenceCCF = 0.05,
  minSnvNum = 5
)
}
\arguments{
\item{ccfInputList}{The 'mutMulti' list in the list returned by the runTumorClutter()
or postProcess() function, both functions return the same 'mutMulti'.}

\item{sampleName}{A vector of sample names, the length of which should be the
same as the number of samples, or the names of the elements of the "ccfInputList"
if not specified.}

\item{minPresenceCCF}{Minimum CCF in the sample to mark SNV presence,
below which SNV will be marked as absent, default 0.05.}

\item{minSnvNum}{Minimum number of SNVs in the group, groups with SNVs less
than this value will be deleted, default is 5, if you don't want to delete
the group you can set it to 0.}
}
\value{
Returns a list containing a grouped data.frame and an average CCF
data.frame for each group.
}
\description{
Group SNVs according to their presence or absence in the sample
and remove groups with too few mutations in the group.
}
