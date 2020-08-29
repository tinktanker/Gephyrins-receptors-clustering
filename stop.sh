#!/bin/bash
stop=1000  #let the code shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c pause every 'stop' steps
tot=2
sh "executei.sh"
for ((n=0;n<tot;n++))   #total=tot*stop
{
  echo "n=$n!!!"
  let A=$n*$stop
  let B=($n+1)*$stop
  sed -i "s/tmax $A/tmax $B/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
  sh "executesg.sh" #run c code
  sed -i "s/tau=$A/tau=$B/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
  sed -i "s/gephyrin$A/gephyrin$B/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
  sed -i "s/receptor$A/receptor$B/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
  sed -i "s/cluster$A/cluster$B/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
  sed -i "s/stat$A/stat$B/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c


  #delete record files which was just usded in this loop
  rm gephyrin$A
  rm receptor$A
  rm cluster$A
  rm stat$A
}
rm gephyrin$B
rm receptor$B
rm cluster$B
rm stat$B

#change the parameters' value in the code back to be the original
sed -i "s/tau=$B/tau=0/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
sed -i "s/tmax $B/tmax 0/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
sed -i "s/gephyrin$B/gephyrin0/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
sed -i "s/receptor$B/receptor0/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
sed -i "s/cluster$B/cluster0/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
sed -i "s/stat$B/stat0/g" shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
