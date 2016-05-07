# Sepsis-Omics Tutorials

Tutorials on how to use the microGVL to do common bacterial bioinformatic tasks.

## Deployment

The tutorials have been deployed here: http://sepsis-omics.github.io/tutorials/

## How to work on them locally

1. Clone the repo
```
% git clone https://github.com/sepsis-omics/tutorials.git sepsis-tutorials
% cd sepsis-tutorials
```

2. Browse the site locally without deploying to public internet
```
% mkdocs serve
```
Open your web browser to http://127.0.0.1:8000/ and leave it open. 
This will update automatically as you make changes to the documenation.

3. The master document is a YAML file
```
% less mkdocs.yml
```

4. The actual Markdown pages are in the `docs` folder:
```
% ls docs
index.md about.md     # some pages
dna/ prot/ rna/ met/  # some folders with more pages
media/                # folder for images
```

5. To add a new page, say a page on the 'Minia' genome assembler, find the right location and create a page. 
In this case it would be `docs/dna/denovo/minia.md`. Write the tutorial in that file, and then add the file to the 
master document `mkdocs.yml` in the correct section.

6. When you are happy, add it the repo:
```
git add docs/dna/denovo/minia.md
git commit -m "Added minia" mkdocs.yml docs/dna/denovo/minia.md
git push
```
Your private local web version  http://127.0.0.1:8000/ will also update.

7. To deploy the whole lot to the *public* website
```
mkdocs gh-deploy --clean --message "Added minia"
```
This first builds a web HTML version of our Markdown hierarchy into the `site/` folder, then pushes it to a special
branch of the github repo called `gh-pages` which GitHub makes available at the public URL http://sepsis-omics.github.io/tutorials/

## Authors

* Torsten Seemann
* Anna Syme
* Madison Flannery
* Simon Gladman
