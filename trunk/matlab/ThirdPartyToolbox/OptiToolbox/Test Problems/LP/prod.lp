\Problem name: prod.lp

Maximize
 obj: - 10 make(1,1) - 10 make(1,2) - 10 make(1,3) - 10 make(1,4)
      - 11 make(2,1) - 11 make(2,2) - 11 make(2,3) - 11 make(2,4)
      + 25 sell(1,1) + 26 sell(1,2) + 27 sell(1,3) + 27 sell(1,4)
      + 30 sell(2,1) + 35 sell(2,2) + 37 sell(2,3) + 39 sell(2,4)
Subject To
 c1:  0.005 make(1,1) + 0.00714285714285714 make(2,1) <= 40
 c2:  0.005 make(1,2) + 0.00714285714285714 make(2,2) <= 40
 c3:  0.005 make(1,3) + 0.00714285714285714 make(2,3) <= 32
 c4:  0.005 make(1,4) + 0.00714285714285714 make(2,4) <= 40
 c5:  make(1,1) - inv(1,1) - sell(1,1)  = -10
 c6:  make(1,2) + inv(1,1) - inv(1,2) - sell(1,2)  = 0
 c7:  make(1,3) + inv(1,2) - inv(1,3) - sell(1,3)  = 0
 c8:  make(1,4) + inv(1,3) - inv(1,4) - sell(1,4)  = 0
 c9:  make(2,1) - inv(2,1) - sell(2,1)  = 0
 c10: make(2,2) + inv(2,1) - inv(2,2) - sell(2,2)  = 0
 c11: make(2,3) + inv(2,2) - inv(2,3) - sell(2,3)  = 0
 c12: make(2,4) + inv(2,3) - inv(2,4) - sell(2,4)  = 0
Bounds
 0 <= sell(1,1) <= 6000
 0 <= sell(1,2) <= 6000
 0 <= sell(1,3) <= 4000
 0 <= sell(1,4) <= 6500
 0 <= sell(2,1) <= 4000
 0 <= sell(2,2) <= 2500
 0 <= sell(2,3) <= 3500
 0 <= sell(2,4) <= 4200
End
