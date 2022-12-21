import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def main(args):
  mm = scipy.io.mmread(args.mm).T.tocsc()
  
  scrub = scr.Scrublet(mm, expected_doublet_rate = args.expDoubletRate)
  doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts = 50, 
    min_cells = 3, 
    min_gene_variability_pctl = args.gene_variation_percentile, 
    n_prin_comps = args.npc)
  
  scrub.call_doublets(threshold = 0.25)
  
  fig, ax = scrub.plot_histogram()
  fig.savefig(args.outDir + "/" + args.sampleName + "_hists.png")
  
  with open(args.outDir + "/" + args.sampleName + "_doubletScores.txt", 'w') as file:
    for line in doublet_scores:
      file.write("{}\n".format(line))


if __name__ == '__main__':
  # set up command line arguments
  parser = argparse.ArgumentParser(
    description = "Identify doublets using scrublet")

  parser.add_argument("--sampleName",
    required = True,
    help = "Sample name")
  parser.add_argument("--mm",
    required = True,
    help = "Raw count matrix in MM format")
  parser.add_argument("--outDir",
    required = True,
    help = "Out directory path")
  parser.add_argument("--expDoubletRate",
    default = 0.06,
    type = float,
    help = "Expected doublet rate. Default is 0.06")
  parser.add_argument("--gene_variation_percentile",
    default = 85,
    type = int,
    help = "Percentile filter for features by variation")
  parser.add_argument("--npc",
    default = 30,
    type = int,
    help = "Number of principal components to use")

  args = parser.parse_args()

  if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)

  if not os.path.exists(args.mm):
    raise Exception("Matrix file not found")


  main(args)
