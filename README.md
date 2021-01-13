# photo_rate_constant_demo
Demonstration of runtime-initialized photodecompostion rate constants in a docker container

Install Docker on your machine
Then git clone this repository and change into the directory
```
docker build -t photo_rate_constant_demo .
docker run -it photo_rate_constant_demo bash
```

Now you are in the container
```
cd photo-demo/
./src/build/photo_rate_demo photo.config.json
```
