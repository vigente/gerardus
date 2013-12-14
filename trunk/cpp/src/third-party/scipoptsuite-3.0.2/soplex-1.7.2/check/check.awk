#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*   File....: check.awk                                                     *
#*   Name....: SoPlex Check Report Generator                                 *
#*   Author..: Thorsten Koch                                                 *
#*   Copyright by Author, All rights reserved                                *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
function abs(x)
{
   if (x < 0)
      return -x;
   if (x > 0)
      return x;
   return 0.0;  # get rid of -0.0
}
function printviol(x)
{
   if (x < 1e-9)
      printf("         ");
   else
      printf(" %.2e", abs(x));
}
BEGIN {
   print "";
   line = "-----------------------------------------------------------------------------------------------------------------------------\n";
   printf(line);
   printf("Name         Rows   Cols Type   Iter     Time Objective        relErr   maxVCons sumVCons maxVBoun sumVBoun maxVRedc sumVRedc\n");

   obj      = "error";
}
/=opt=/          { sol[$2] = $3; }
/=type=/         { type = $2; }
/IEXAMP22/       { file = $5; }
/IEXAMP24/       { rows = $4; cols = $6; }
/Solution time/  { time = $4; } 
/Iterations/     { iter = $3; }
/IEXAMP29/       { obj  = $5; }
/IEXAMP31/       { infeas = 1; }
/IEXAMP32/       { infeas = 1; } 
/IEXAMP33/       { timeout = 1; }
#/EEXAMP39/       { singular = 1; }
/XSOLVE21/       { singular = 1; }
#/EEXAMP40/       { cycling = 1; }
/XSOLVE13/       { cycling = 1; }
/IEXAMP07/       { cvm = $4; cvs = $5; if (cvm > cvmax[type]) cvmax[type] = cvm; cvsum[type] += cvs; }
/IEXAMP09/       { bvm = $4; bvs = $5; if (bvm > bvmax[type]) bvmax[type] = bvm; bvsum[type] += bvs; }
/IEXAMP11/       { rcm = $5; rcs = $6; if (rcm > rcmax[type]) rcmax[type] = rcm; rcsum[type] += rcs; }
/=start=/        {
   type = "";
   for(i = 2; i <= NF; i++)
      type = type substr($i, 2);
}
/ready/       {
   n = split(file, a, "/");
   split(a[n], b, ".");
   name = b[1];

   if (sol[name] == "")
      print name, "nicht gefunden";
   else
   {
      if (name == prevname)
	 printf("%25s", "");
      else
      {
	 printf(line);
	 printf("%-10s %6d %6d ", name, rows, cols);
      }
      printf("%-3s %7d %8.2f ", type, iter, time);

      if (infeas)
	 printf("%-14s", "infeasible");
      else if (timeout)
	 printf("%-14s", "timeout");
      else if (cycling)
	 printf("%-14s", "cycling");
      else if (singular)
	 printf("%-14s", "singular");
      else if (obj == "error")
	 printf("%-14s", "error");
      else
	 printf("%+e ", obj);

      if (timeout)
	 printf("\n");
      else if ( obj == "error" && !infeas)
      {
	 printf("XX\n");
	 fail[type]++;
	 fails++;
      }
      else
      {
	 if (!infeas && sol[name] != "infeasible")
	 {
	    abserr = abs(sol[name] - obj);
	    if (abs(sol[name]) >= 1e-5)
	       relerr = abserr / abs(sol[name]);
	    else
	       relerr = abserr;

	    if ((abserr < 1e-4) || (relerr < 1e-5))
	    {
	       printf("ok %.2e", relerr);
	       pass[type]++;
	       relerrsum[type] += relerr;
	       passes++;
	    }
	    else
	    {
	       printf("XX %.2e", abserr);
	       fail[type]++;
	       fails++;
	    }
	    printviol(cvm);
	    printviol(cvs);
	    printviol(bvm);
	    printviol(bvs);
	    printviol(rcm);
	    printviol(rcs); 
	    print "";
	 }
	 else
	 {
	    if (infeas == 1 && sol[name] == "infeasible")
	    {
	       printf("ok\n");
	       pass[type]++;
	       passes++;
	    }
	    else
	    {
	       if (infeas && sol[name] != "infeasible")
		  printf("XX %.2e\n", abs(sol[name]));
	       else
		  printf("XX infeasible\n");
		  
	       fail[type]++;
	       fails++;
	    }
	 }
      }
      sum[type] += time;
      cnt[type]++;
      counts++;
      times += time;
   }
   prevname = name;
   timeout  = 0;
   infeas   = 0;
   singular = 0;
   cycling  = 0;
   obj      = "error";
   iter     = 0;
   time     = 0;
   rows     = 0;
   cols     = 0;
}
END {
   print "";
   printf(line);
   printf("Alg            Cnt  Pass  Fail       Time                      relErr   maxVCons sumVCons maxVBoun sumVBoun maxVRedc sumVRedc\n");
   printf(line);
   for(i in sum)
   {
      printf("%-12s %5d %5d %5d %10.2f                     ",
	     i, cnt[i], pass[i], fail[i], sum[i]);
      printviol(relerrsum[i]);
      printviol(cvmax[i]);
      printviol(cvsum[i]);
      printviol(bvmax[i]);
      printviol(bvsum[i]);
      printviol(rcmax[i]);
      printviol(rcsum[i]);
      print "";

      relerrsumsum += relerrsum[i];
      cvmaxsum += cvmax[i];
      cvsumsum += cvsum[i];
      bvmaxsum += bvmax[i];
      bvsumsum += bvsum[i];
      rcmaxsum += rcmax[i];
      rcsumsum += rcsum[i];
   }
   printf(line);
   printf("%-12s %5d %5d %5d %10.2f                     ",
	  "Sum", counts, passes, fails, times);

   printviol(relerrsumsum);
   printviol(cvmaxsum);
   printviol(cvsumsum);
   printviol(bvmaxsum);
   printviol(bvsumsum);
   printviol(rcmaxsum);
   printviol(rcsumsum);
   print "";
   printf(line);
}






