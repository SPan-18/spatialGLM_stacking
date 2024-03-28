rm(list = ls())

# merge BBS data 2010-2019

bbs10 <- read.csv("bbs/BBS10.csv")
bbs11 <- read.csv("bbs/BBS11.csv")
bbs12 <- read.csv("bbs/BBS12.csv")
bbs13 <- read.csv("bbs/BBS13.csv")
bbs14 <- read.csv("bbs/BBS14.csv")
bbs15 <- read.csv("bbs/BBS15.csv")
bbs16 <- read.csv("bbs/BBS16.csv")
bbs17 <- read.csv("bbs/BBS17.csv")
bbs18 <- read.csv("bbs/BBS18.csv")
bbs19 <- read.csv("bbs/BBS19.csv")
# bbs21 <- read.csv("BBS21.csv")

# bbs21 <- bbs21 %>% 
#   filter(RouteDataID == 6376050)

bbs10$Time <- 1
bbs11$Time <- 2
bbs12$Time <- 3
bbs13$Time <- 4
bbs14$Time <- 5
bbs15$Time <- 6
bbs16$Time <- 7
bbs17$Time <- 8
bbs18$Time <- 9
bbs19$Time <- 10
# bbs21$Time <- 6



bbs <- do.call("rbind", list(bbs10, bbs11, bbs12, bbs13, bbs14, 
                             bbs15, bbs16, bbs17, bbs18, bbs19))
bbs <- bbs[, -1]

write.csv(bbs, "BBS10-19.csv", row.names = F)
