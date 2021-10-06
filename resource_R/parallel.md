```
# load packages
if(!require(parallel)) {
  install.packages("parallel")
}
library(parallel)

# 코어 개수 획득
numCores <- parallel::detectCores() - 1

# 클러스터 초기화
myCluster <- parallel::makeCluster(numCores)

# 변수 등록 (important; to share the variable with multiple cores)
base <- 2
parallel::clusterExport(myCluster, "base")

# CPU 병렬처리
parallel::parLapply(cl = myCluster,
                    X = 2:4,
                    fun = function(x) {
                      base^x
                    })

# 클러스터 중지
parallel::stopCluster(myCluster)
```

# parSapply
```
# load packages
if(!require(parallel)) {
  install.packages("parallel")
}
library(parallel)

# 코어 개수 획득
numCores <- parallel::detectCores() - 1

# 클러스터 초기화
myCluster <- parallel::makeCluster(numCores, type = "PSOCK")
myCluster <- parallel::makeCluster(numCores, type = "FORK")

# 변수 등록
base <- 2
parallel::clusterExport(myCluster, "base")

# CPU 병렬처리
parSapply(myCluster, as.character(2:4), 
          function(exponent){
            x <- as.numeric(exponent)
            c(base = base^x, self = x^x)
          })

# 클러스터 중지
parallel::stopCluster(myCluster)
```



# foreach
foreach 패키지는 루프(Loop)와 lapply() 함수를 융합한 것으로 R의 병렬처리에 있어 매우 인기가 많다.
foreach() 함수의 .combine 속성은 결과를 어떻게 결합시킬 것인가를 정의한다.
가령, Vector 형태의 결과를 출력하려면

```
# load packages
if(!require(foreach)) {
  install.packages("foreach")
}
library(foreach)

if(!require(doParallel)) {
  install.packages("doParallel")
}
library(doParallel)

# 코어 개수 획득
numCores <- parallel::detectCores() - 1

# 클러스터 초기화
myCluster <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(myCluster)

# 변수 등록
base <- 2
parallel::clusterExport(myCluster, "base")

# CPU 병렬처리
foreach::foreach(exponent = 2:4, .combine = c)  %dopar% {
  base^exponent
}

# 클러스터 중지
parallel::stopCluster(myCluster)
```
