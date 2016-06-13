# Phandango

FIXME: current Firefox on mGVL Lubuntu desktop is too old to use Phandango


FIXME: ideally: sepsis related: would be good to see clade x has AMR genes? e.g. this example shows a subclade which is ceftriaxone susceptable but azithromycin resistant https://github.com/jameshadfield/phandango/wiki/Panel%20Layout#metadata

Phandango is a tool to vizualize genome phylogenies.


## Start

- prepare input files (e.g. from Roary output):
    - <fn>gene_presence_absence.csv</fn>: gene name, various information, then a column for each sample - if the gene is present in that sample, it is listed here with an appended gene ID.
    - <fn>tree.newick</fn>: a phylogenetic tree based on an alignment of core genes. (Note: this is not a default Roary output and must be specified by creating an alignment of core genes, and then building a phylogenetic tree).

- navigate to these files in your mGVL (e.g. ssh in terminal) and move them to your public_html folder.
- open your public_html folder (e.g. http://mgvl_IP/public/username  or similar), and download to your local computer.

## Run

- open the Phandango [webpage.](http://jameshadfield.github.io/phandango/)

- drag and drop these files onto the webpage.

- now we can see the tree on the left, and the core and accessory genome aligned to each sample.

![Phandango screenshot](./images/image00.png)

- things to look at:

    - top tabs: settings - change the displayed labels
    - change panel sizes - drag grey circles at the edges of each panel
    - line graph?


## Save output image

- press p to save the displayed data as a vector SVG file.


## What next
