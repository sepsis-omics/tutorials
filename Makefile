
.PHONY: all strict serve deploy

all:
	mkdocs build --clean

strict:
	mkdocs build --clean --strict

serve:
	mkdocs serve

deploy:
	mkdocs gh-deploy --clean --message "publicly deploy"
	