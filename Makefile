
all: ui/guiView.py ui/InputTab.py

ui/guiView.py: ressources/guiView.ui
	pyuic4 ressources/guiView.ui > $@

ui/InputTab.py: ressources/InputTab.ui
	pyuic4 ressources/InputTab.ui > $@
