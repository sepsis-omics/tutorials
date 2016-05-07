
.PHONY: all strict serve

all:
	mkdocs build --clean

strict:
	mkdocs build --clean --strict

serve:
	mkdocs serve
		