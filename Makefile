docker:
	docker run -p 8888:8888 --rm --mount type=bind,source=$$(pwd),target=/home/jovyan/work fermi-notebook

build:
	docker build -t fermi-notebook .

dot:
	gprof2dot -f pstats fermi.pstats | dot -Tpng -o fermi_profile.png
