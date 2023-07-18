include .entangled/makefile.in

.PHONY: figures deploy-pages

figures: figures.mk
> make -f figures.mk

deploy-pages:
> git checkout gh-pages; \
> rsync -r site/* .; \
> git commit -a -m 'deploy pages'; \
> git push origin gh-pages

