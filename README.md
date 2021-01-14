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
 photo_rate_demo: Photorate constants at the top of the atmosphere (1/s)
  2.1640320E-05  5.9827795E-06  5.6364924E-05  3.8730194E-04  2.1947110E-03
  3.2114790E-04  1.4738591E-04  2.5643905E-03
```
