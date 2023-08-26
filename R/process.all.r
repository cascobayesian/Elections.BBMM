kl.dist <- function(x,y)
{
	if (!is.na(x)&!is.na(y))
	{
		if ((x > 0)&(y>0))
		{	
			return(x*log(x/y))
		}
		else
		{
			return(NA)
		}
	}else
	{
		return(NA)
	}
}

js.dist <- function(x,y)
{
	return(0.5*kl.dist(x,y)+0.5*kl.dist(y,x))
}

calc.confusion <- function(clusters)
{
  T<- dim(clusters)[1]
  N <- dim(clusters)[2]

  confusion <- array(0,dim=c(N,N))

  for(i in 1:(N-1))
  {
  	for(j in (i+1):(N))
  	{
    	confusion[i,j] <- sum( clusters[,i]==clusters[,j] )/T;
    	confusion[j,i] <- confusion[i,j]
  	}
  }
  return(confusion)
}

calc.probs <- function(confusion,r)
{
	K <- max(r)
	p <- array(0,dim=c(K,length(r)))
	o <- order(r)

	for(s in 1:length(r))
	{
		tmp <- confusion[o,o][s,]
		cs <- cumsum(table(r))
		
		for(t in 1:K)
		{
			if (t==1)
			{
				p[t,s] <- quantile(tmp[1:cs[1]])[3]
			}else
			{
				p[t,s]<- quantile(tmp[(cs[t-1]+1):cs[t]])[3] 
			}
		}
		p[,s] <- p[,s]/sum(p[,s])
	}
	return(p)
}

draw.circle <- function(x,center,r=0.5,cols)
{
  K <- length(x)
  
  cirlce <- seq(0,2*pi,by=0.01)
  #cols <- rainbow(K)
  cs <- c(0,cumsum(x));
  for(i in 1:(K))
  {
    x <- c(0,cos(seq(2*pi*cs[i],2*pi*cs[i+1],by=0.001)),0)*r+center[1]
    y <- c(0,sin(seq(2*pi*cs[i],2*pi*cs[i+1],by=0.001)),0)*r+center[2]
    
    polygon(x,y,col=cols[i],border=cols[i])
  }
}

## processing


## set path/prefix

## load 
load('Output/Multnomah.Ref/multnomah_later_processed_3.RData')
#load('Output/Maine.Ref/out.1.all.RData')
## pull out clusters and times
results.file <- "test.csv"#"mass_k=6_run=2.csv"
figure.file <-  'fig.pdf'#"mass_run=2_k=6.pdf"
pdf(figure.file,width=12,height=9)
N <- dim(dataSet)[1]
Q <- dim(dataSet)[2]

cat("N:",N,"\t Q:",Q,"\n")

clusters <- array(0,dim=c(T/2,N))
times <- array(0,dim=c(T/2,3))

base <- T/2 - 1
for(t in 1:(T/2))
{
	r <- outs[[t+base]][[1]]$r
	clusters[t,] <- r
	times[t,] <- c(outs[[t+base]][[2]],sum(table(r)>=1),sum(table(r)>=0))
}

## calcuate distribution of K
max.K <- max(times[,2])
min.K <- min(times[,2])
seq.K <- min.K:max.K

k.w <- rep(0,length(seq.K))
names(k.w) <- seq.K
for(k in 1:length(seq.K))
{
	k.pos <- which(times[,2]==seq.K[k])

	k.w[k] <- sum(times[k.pos,1])/sum(times[,1])
}

plot(seq.K,k.w,type='h',lwd=3,axes=FALSE,xlab="Number of clusters",ylab="Posterior Probability")
axis(1,at=seq.K,labels=seq.K)
axis(2)
box()

## calc posterior mode K

post.k<- seq.K[which.max(k.w)]

## calc confusion matrix

confusion.w <- array(0,dim=c(N,N))
confusion.u <- array(0,dim=c(N,N))

for(i in 1:(N-1))
{
	confusion.w[i,i] <- 1
	confusion.u[i,i] <- 1
	for(j in (i+1):N)
	{
		confusion.w[i,j] <- sum((clusters[1:(T/2),i]==clusters[1:(T/2),j])*times[1:(T/2),1])/sum(times[1:(T/2),1])
		confusion.w[j,i] <- confusion.w[i,j]
		confusion.u[i,j] <- sum(clusters[1:(T/2),i]==clusters[1:(T/2),j])/(T/2)
		confusion.u[j,i] <- confusion.u[i,j]
	}

}
confusion.w <- 2*confusion.w-1
confusion.u <- 2*confusion.u-1

## calc representative clustering

library(cluster)
c.w <- pam(confusion.w,post.k)$clustering
c.u <- pam(confusion.u,post.k)$clustering

o.w <- order(c.w)
o.u <- order(c.u)

## confusion matrix ordered, with boudaries
cols <- rainbow(100)
par(mgp=c(2,1,0),cex.lab=1.5)
image(confusion.w[o.w,o.w],axes=FALSE,xlab='Community (by voting bloc)',ylab='Community (by voting bloc)')
box()

usr <- par("usr")
off.set <- usr[1]
scale <- abs(usr[1])+usr[2]
b <- cumsum(table(c.w))/length(c.w)
for(k in 1:(post.k-1))
{
	pos <- (scale*b[k]+off.set)
	abline(v=pos,col=grey(0.9,0.8),lwd=3)
	abline(h=pos,col=grey(0.9,0.8),lwd=3)
}
b.new <- c(0,b)
diff <- (b.new[1:post.k]+b.new[2:(post.k+1)])/2
axis(1,at=diff,labels=1:post.k,col=grey(0.5))
axis(2,at=diff,labels=1:post.k,col=grey(0.5))

## calc barplot

p <- dataSet[,,1]/(dataSet[,,1]+dataSet[,,2])

out.mat <- matrix(0,nrow=1,ncol=2)
out.mat[1,] <- c(1,2)
layout(out.mat,heights=rep(3,2),widths=c(3,1))
image(p[o.w,],col=cols,xlab='Community (by voting bloc)',ylab="Questions (by year)",axes=FALSE)
box()
for(k in 1:(post.k-1))
{
	pos <- (scale*b[k]+off.set)
	abline(v=pos,col=grey(0.9,0.8),lwd=3)
}
plot.new()
q <- round(quantile(p,c(0.05,0.25,0.5,0.75,0.95),na.rm=TRUE),2)
legend(x='center',legend=q,col=cols[c(5,25,50,75,95)+1],pch=20,bg='white',cex=1.5,title='Support')

p.new <- p
for(q in 1:Q)
{
	p.new[,q] <- (p.new[,q]-mean(p.new[,q],na.rm=TRUE))/sd(p.new[,q],na.rm=TRUE)
}

out.mat <- matrix(0,nrow=1,ncol=2)
out.mat[1,] <- c(1,2)
layout(out.mat,heights=rep(3,2),widths=c(3,1))
image(p.new[o.w,],col=cols,xlab='Community (by voting bloc)',ylab="Questions (by year)",axes=FALSE)
box()
for(k in 1:(post.k-1))
{
	pos <- (scale*b[k]+off.set)
	abline(v=pos,col=grey(0.9,0.8),lwd=3)
}

plot.new()
q <- round(quantile(p.new,c(0.05,0.25,0.5,0.75,0.95),na.rm=TRUE),2)
legend(x='center',legend=q,col=cols[c(5,25,50,75,95)+1],pch=20,bg='white',cex=1.75,title='Z-score')

## calc probs
probs <- calc.probs(0.5*(confusion.u+1),c.u)
image(log(probs),col=cols,axes=FALSE,xlab='Voting bloc')
#axis(1)
box()
for(k in 1:(post.k-1))
{
	pos <- (scale*(k/(post.k-1))+off.set)
	#abline(v=pos,col=grey(0.9,0.8),lwd=3)
}
plot.new()
q <- round(quantile(log(probs),c(0.05,0.25,0.5,0.75,0.95)),2)
legend(x='center',legend=round(exp(q),2),col=cols[c(5,25,50,75,95)+1],pch=20,bg='white',cex=1.75,title='Probability')

print('hi')
#dev.off()
## calc tsne

library(Rtsne)
## support by question	

#bad.town <- c(8)
#N <- N-length(bad.town)

#c.w <- c.w[-c(bad.town)]
#o.w <- order(c.w)

#towns.un <- towns[-bad.town]

#dataSet <- dataSet[-bad.town,,]

js.all <- array(0,dim=c(N,N,Q))
js.mat <- array(0,dim=c(N,N))
a.all <- js.all*0;
a.mat <- js.mat*0;

for(j in 1:(N-1))
{
	

	p <- dataSet[j,,1]/(dataSet[j,,1]+dataSet[j,,2])
	y <- log(dataSet[j,,1]/dataSet[j,,2])
	print(j)
	for(i in (j+1):N)
	{
		x <- log(dataSet[i,,1]/dataSet[i,,2])
		q <- dataSet[i,,1]/(dataSet[i,,1]+dataSet[i,,2])		

		for(k in 1:Q)
		{	
			
			if ((x[k]!=Inf)&(y[k]!=Inf)&(!is.na(x[k]))&(!is.na(y[k])))
			{
				a.all[i,j,k] <- dist(c(x[k],y[k]))
				a.all[j,i,k] <- a.all[i,j,k]
				
			}else
			{
				a.all[i,j,k] <- NA
				a.all[j,i,k] <- NA

			}
			js.all[i,j,k] <- js.dist(p[k],q[k])
			js.all[j,i,k] <- js.all[i,j,k]
		}
		a.mat[i,j] <- mean(a.all[i,j,],na.rm=TRUE)
		a.mat[j,i] <- mean(a.all[j,i,],na.rm=TRUE)
		js.mat[i,j] <- mean(js.all[i,j,],na.rm=TRUE)
		js.mat[j,i] <- js.mat[i,j]
	}
}

print('hi')
cond.prob <- array(0,dim=c(N,N))

s.2 <- 0.05
for(j in 1:(N-1))
{

	for(i in (j+1):N)
	{
		cond.prob[j,i] <- exp(-js.mat[j,i]/(2*s.2))/sum(exp(-js.mat[j,-c(i)]/2*s.2),na.rm=TRUE)
		cond.prob[i,j] <- cond.prob[j,i]
	}
}


test <- Rtsne(cond.prob,dims=2,perplexity=10)

conf <- (confusion.w+1)/2
probs <- calc.probs(conf,c.w)
par(cex.lab=1.5,mar=c(5.5,5,4,2),mfrow=c(1,2))
clust <- apply(probs,2,which.max)

cols <- rainbow(post.k,s=0.5,alpha=0.75)
plot(test$Y[o.w,1],test$Y[o.w,2],pch=20,xlab='Component 1', ylab = 'Component 2',col=cols[clust],cex=2)
plot(test$Y[o.w,1],test$Y[o.w,2],pch=20,xlab='Component 1', ylab = 'Component 2',col="white")

for(j in c(1:N))
{
	center <- c(test$Y[o.w[j],1],test$Y[o.w[j],2])
	draw.circle(probs[,j],center,r=0.75,cols)
	text(center[1],center[2],labels=towns[o.w[j]],cex=0.25,col=grey(0.5,0.5))

}
legend(x='topright',legend=1:post.k,pch=20,col=cols,cex=1.5)

dev.off()

overall.support <- array(0,dim=c(Q,2*post.k))

support <- t(dataSet[,,1]/(dataSet[,,1]+dataSet[,,2]))

odd <- seq(1,2*post.k,by=2)
even <- seq(2,2*post.k,by=2)

for(q in 1:Q)
{
	for (k in 1:post.k)
	{
		overall.support[q,2*k-1] <- mean(support[q,c.w==k])
		overall.support[q,2*k] <- sd(support[q,c.w==k])
	}
}

col.names <- rep(c("Mean","SD"),post.k)
for(k in 1:post.k)
{
	col.names[2*k-1] <- paste(col.names[2*k-1],"_bloc_",k,sep="")
	col.names[2*k] <- paste(col.names[2*k],"_bloc_",k,sep="")
}
row.names <- paste("Q",1:Q,sep="")

colnames(overall.support) <- col.names
rownames(overall.support) <- row.names

write.table(file=results.file,x=round(overall.support,3),col.names=TRUE,row.names=TRUE,sep=",")