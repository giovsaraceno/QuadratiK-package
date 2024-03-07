library(testthat)
test_that("Clustering algorithm works", {
   
   #------------------------------------------------------
   ## Clustering on the Sphere
   ## 
   size <- 100
   groups<-c(rep(1, size), rep(2, size),rep(3,size))
   rho=0.8
   data1<-rpkb(size, c(1,0,0),rho,method='rejvmf')
   data2<-rpkb(size, c(0,1,0),rho,method='rejacg')
   data3<-rpkb(size, c(-1,0,0),rho,method='rejpsaw')
   dat<-rbind(data1$x,data2$x, data3$x)
   expect_equal(dim(dat),c(size*3,3))
   
   pkbd_res<- pkbc(dat, 3)
   
   expect_true(class(pkbd_res)== "pkbc")
   expect_type(pkbd_res@res_k, "list")

})
