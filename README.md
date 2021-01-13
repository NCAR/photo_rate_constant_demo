# photo_rate_constant_demo
Demonstration of runtime-initialized photodecompostion rate constants in a docker container

Install Docker on your machine

Clone this repository and change into the resulting directory and execute the following
```
docker build -t photo_rate_constant_demo .
docker run -it photo_rate_constant_demo bash
```

Now you are in the container.
```
cd photo-demo/
./src/build/photo_rate_demo photo.config.json
```

The expected output follows
```
 photo_rate_demo: Opened and read config file photo.config.json
  
 Photo kinetics constructor: key = CH3CHO+hv->CH3+HCO
 cross section builder: entering
 base cross_section constructor: entering
 ERROR:initialize /Users/stacy/Documents/Python_dev/Sandbox/XSQY/Data/XSQY/CH3CHO_cross_section.nc could not be created
 read_netcdf_file: retcode from initialize =            2
STOP FileOpenError
```
