# Writeup

I found that using a k of 3 and a p value of 2 was the best option.  When I increased the value of k I recieved a lower accuracy. This decrease in accuracy was almost linear. When I tried a k value of 1, I also recieved a lower accuracy. I tried many values, up to k = 111. Therefore, while it feels wrong to say the default value, from my results, a k value of 3 was best. In terms of p values, I found that using a p value of 2 resulted in better accuracy than using a p value of 2. While using a p value of 3 gave a slightly more accurate result, it took much much longer and I don't think that the slight increase in accuracy was worth the increase in time. 


# Results

## k = 1, p = 2
For a k-NN classifier using 1 neighbors and a Minkowski distance of order 2
Confusion Matrix (Truth \ Predicted):
           0    1    2    3    4    5    6    7    8    9
       --------------------------------------------------
    0 |  973    1    1    0    0    1    3    1    0    0
    1 |    0 1129    3    0    1    1    1    0    0    0
    2 |    7    6  992    5    1    0    2   16    3    0
    3 |    0    1    2  970    1   19    0    7    7    3
    4 |    0    7    0    0  944    0    3    5    1   22
    5 |    1    1    0   12    2  860    5    1    6    4
    6 |    4    2    0    0    3    5  944    0    0    0
    7 |    0   14    6    2    4    0    0  992    0   10
    8 |    6    1    3   14    5   13    3    4  920    5
    9 |    2    5    1    6   10    5    1   11    1  967
Classification accuracy = 96.91%

## k = 3, p = 2
For a k-NN classifier using 3 neighbors and a Minkowski distance of order 2
Confusion Matrix (Truth \ Predicted):
           0    1    2    3    4    5    6    7    8    9
       --------------------------------------------------
    0 |  974    1    1    0    0    1    2    1    0    0
    1 |    0 1133    2    0    0    0    0    0    0    0
    2 |   10    9  996    2    0    0    0   13    2    0
    3 |    0    2    4  976    1   13    1    7    3    3
    4 |    1    6    0    0  950    0    4    2    0   19
    5 |    6    1    0   11    2  859    5    1    3    4
    6 |    5    3    0    0    3    3  944    0    0    0
    7 |    0   21    5    0    1    0    0  991    0   10
    8 |    8    2    4   16    8   11    3    4  914    4
    9 |    4    5    2    8    9    2    1    8    2  968
Classification accuracy = 97.05%

## k = 3, p = 1
For a k-NN classifier using 3 neighbors and a Minkowski distance of order 1
Confusion Matrix (Truth \ Predicted):
           0    1    2    3    4    5    6    7    8    9
       --------------------------------------------------
    0 |  974    1    1    0    0    1    2    1    0    0
    1 |    0 1134    1    0    0    0    0    0    0    0
    2 |   11   16  979    3    1    0    1   17    4    0
    3 |    1    3    4  974    1   14    0    9    2    2
    4 |    1   13    0    0  934    0    6    4    0   24
    5 |    5    2    0   10    3  859    5    1    1    6
    6 |    6    3    0    0    4    3  942    0    0    0
    7 |    0   25    4    0    2    0    0  988    0    9
    8 |   11    4    6   20    8   17    4    5  895    4
    9 |    3    7    2    9   10    3    1   18    2  954
Classification accuracy = 96.33%

## k = 3, p = 3
For a k-NN classifier using 3 neighbors and a Minkowski distance of order 3
Confusion Matrix (Truth \ Predicted):
           0    1    2    3    4    5    6    7    8    9
       --------------------------------------------------
    0 |  975    1    1    0    0    1    1    1    0    0
    1 |    0 1133    2    0    0    0    0    0    0    0
    2 |    9    7 1000    2    0    0    0   12    2    0
    3 |    0    0    4  980    1   11    1    7    3    3
    4 |    2    6    0    0  948    0    4    2    0   20
    5 |    4    1    0    9    2  861    6    1    3    5
    6 |    4    3    0    0    3    3  945    0    0    0
    7 |    0   18    6    1    2    0    0  990    0   11
    8 |    6    2    4   17    7    8    3    4  919    4
    9 |    4    6    2    8    9    2    1    9    1  967
Classification accuracy = 97.18%

## k = 5, p = 2
For a k-NN classifier using 5 neighbors and a Minkowski distance of order 2
Confusion Matrix (Truth \ Predicted):
           0    1    2    3    4    5    6    7    8    9
       --------------------------------------------------
    0 |  974    1    1    0    0    1    2    1    0    0
    1 |    0 1133    2    0    0    0    0    0    0    0
    2 |   11    8  991    2    1    0    1   15    3    0
    3 |    0    3    3  976    1   13    1    6    3    4
    4 |    3    7    0    0  944    0    4    2    1   21
    5 |    5    0    0   12    2  862    4    1    2    4
    6 |    5    3    0    0    3    2  945    0    0    0
    7 |    0   22    4    0    3    0    0  988    0   11
    8 |    8    3    5   13    6   12    5    5  913    4
    9 |    5    7    3    9    7    3    1   10    2  962
Classification accuracy = 96.88%

## k = 10, p = 2
For a k-NN classifier using 10 neighbors and a Minkowski distance of order 2
Confusion Matrix (Truth \ Predicted):
           0    1    2    3    4    5    6    7    8    9
       --------------------------------------------------
    0 |  972    1    1    0    0    2    3    1    0    0
    1 |    0 1132    2    0    0    0    1    0    0    0
    2 |   13   12  982    2    1    0    2   17    3    0
    3 |    0    3    3  976    1   10    1    7    6    3
    4 |    2   11    0    0  940    0    4    1    1   23
    5 |    4    0    0   12    1  863    6    1    1    4
    6 |    6    4    0    0    3    2  943    0    0    0
    7 |    0   27    4    0    2    0    0  983    0   12
    8 |    6    4    5   11    7    9    4    7  914    7
    9 |    7    6    3    7   10    3    1   10    2  960
Classification accuracy = 96.65%

## k = 15, p = 2
For a k-NN classifier using 15 neighbors and a Minkowski distance of order 2
Confusion Matrix (Truth \ Predicted):
           0    1    2    3    4    5    6    7    8    9
       --------------------------------------------------
    0 |  970    1    1    0    0    2    5    1    0    0
    1 |    0 1131    2    1    0    0    1    0    0    0
    2 |   15   15  968    3    1    0    3   20    7    0
    3 |    0    3    2  975    1   14    0    7    4    4
    4 |    1   13    0    0  934    0    5    2    1   26
    5 |    3    1    0   10    1  863    8    2    0    4
    6 |    7    4    0    0    3    1  943    0    0    0
    7 |    0   28    3    0    2    0    0  980    0   15
    8 |    7    4    5   13    7   12    5    7  907    7
    9 |    6    7    2    9    9    2    1   10    1  962
Classification accuracy = 96.33%

## k = 45, p = 2
For a k-NN classifier using 45 neighbors and a Minkowski distance of order 2
Confusion Matrix (Truth \ Predicted):
           0    1    2    3    4    5    6    7    8    9
       --------------------------------------------------
    0 |  968    1    1    0    0    2    7    1    0    0
    1 |    0 1130    2    1    0    0    2    0    0    0
    2 |   23   28  942    6    2    0    4   20    7    0
    3 |    0    4    2  972    1   11    1    9    6    4
    4 |    0   16    0    0  927    0    8    1    1   29
    5 |    3    7    0   10    1  850   12    2    0    7
    6 |    7    5    0    0    3    2  941    0    0    0
    7 |    0   35    4    0    2    0    0  968    0   19
    8 |    8    6    4   18    9   13    6    8  894    8
    9 |    8    7    3   10    7    2    1   12    0  959
Classification accuracy = 95.51%

## k = 111, p = 2
For a k-NN classifier using 111 neighbors and a Minkowski distance of order 2
Confusion Matrix (Truth \ Predicted):
           0    1    2    3    4    5    6    7    8    9
       --------------------------------------------------
    0 |  967    1    1    0    0    3    7    1    0    0
    1 |    0 1130    2    1    0    0    2    0    0    0
    2 |   23   47  906    8    2    1    6   30    9    0
    3 |    0    8    3  961    1   10    1   13    7    6
    4 |    0   21    0    0  901    0   11    2    2   45
    5 |    5   10    0   13    1  836   14    2    0   11
    6 |    9    9    0    0    2    3  935    0    0    0
    7 |    0   44    2    0    3    0    0  954    0   25
    8 |   12   12    2   21   12   17    5   10  867   16
    9 |    8   10    2    7    6    3    1   14    0  958
Classification accuracy = 94.15%
