/* 
 * Copyright © 2009 University of Oxford
 *
 * University of Oxford means the Chancellor, Masters and Scholars of
 * the University of Oxford, having an administrative office at
 * Wellington Square, Oxford OX1 2JD, UK. 
 *
 * This file is part of Gerardus.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. The offer of this
 * program under the terms of the License is subject to the License
 * being interpreted in accordance with English Law and subject to any
 * action against the University of Oxford being under the jurisdiction
 * of the English Courts.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */


package gvdecoder;
import java.io.*;
import java.util.*;


import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.geom.*;

import java.awt.image.*;
import java.awt.image.renderable.*;
import javax.media.jai.*;
import javax.media.jai.iterator.*;
import javax.media.jai.widget.*;

import java.awt.color.*;

import java.nio.*;

import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import gvdecoder.log.*;
import gvdecoder.utilities.*;




public class HamamatsuDecoder extends ImageDecoderAdapter {
RandomAccessFile Picfile;
public int FRAMESTEP=1;
int FrameNumber=0;
ArrayList images;
String filename;
String imagepath="";
long frameposition=0;
int orientation=2;
byte[] internalframe=new byte[528];
public static String navfilename;

public int width;
public int height;
public int numframes;
public  int pixelsize;

public boolean baselinesubtract=false;
public int firstindex=0;
public int lastindex=2028;

public  long header_size=50;

public int verMajor;
public int verMinor;
public int sTimeMonth;
public int sTimeDate;
public int sTimeYear;
public int sTimeHour;
public int sTimeMin;
public int sTimeSec;
public int sDataType;
public int sDataByte;
public int sNumChns;
public int sChnStcVer;
public int iNumFrames;
public int sAcqDevice;
public int sChnInt;
public float fScanInt;
public float fCamp;
public float fAmp;
public int sIsAnalyzed;
public int sAnalType;
public int iCommentSize;
public char[] rest=new char[256];

public MappedByteBuffer buffer;
public ShortBuffer shortbuffer;
public short[] singleframe;

public boolean SHOW_V_ONLY=false;
public boolean PAD=false;

public HamamatsuDecoder(){
 FrameNumber=0;
 this.filename=filename;
}


public int CloseImageFile(int instance){return -1;}


public int NumberOfFrames=2028;

//public int FilterOperation(int a, int b, int c, int d){ return 0;}

public int OpenImageFile(String filename){
 int j=0;
 //open the file...
 try{

      Picfile=new RandomAccessFile(filename,"r");
      FileChannel fc = Picfile.getChannel();
      buffer = fc.map (FileChannel.MapMode.READ_ONLY, 0, fc.size());
      buffer.order(ByteOrder.LITTLE_ENDIAN);
      buffer.position(32+1024+4*10+6);
      NumberOfFrames=(buffer.getInt()/512);
      buffer.position(32+1024+140+512);
      shortbuffer=buffer.asShortBuffer();
      singleframe=new short[256];

      return 1;

  }
  catch(EOFException e){System.out.println("eof: number of frames read is "+j);}
  catch(IOException e){System.out.println("error reading the file, please check..."); e.printStackTrace();}
   return 0;

  }

/*
public void stich(){
	//keep header of the first one, load in successive frames

}
*/

public void test(){
	try{
	for (int i=1;i<50;i++){
	Picfile.seek (header_size+256+((520*2)*i));
	System.out.println(getHamamatsuShort(Picfile)+" "+getHamamatsuShort(Picfile));
	}
	}catch(Exception e){e.printStackTrace();}
	}


public int FilterOperation(int OperationCode, int startx, int endx, int instance){
	   if (OperationCode==0){
		   FRAMESTEP=startx;
		   if (FRAMESTEP<1) FRAMESTEP=1;
	   }else
	   if (OperationCode==1){
		   setupbaselinesubtract(startx,endx);
	   }
	   return 0;
	 }
   public String ReturnSupportedFilters(){return("FRAMESTEP,BASELINESUBTRACT");}
   public int ReturnFrameNumber(){return FrameNumber;}
   public int JumpToFrame(int framenum, int instance){
	   //seek to the right location in the file...
       FrameNumber=framenum;
	   shortbuffer.position(FrameNumber*256);

	   return FrameNumber;
	   }


public int ReturnXYBandsFrames(int[] arr, int instance){
	//String fileName=imagepath+File.separator+(String)images.get(FrameNumber);

	//int frames=images.size();

	arr[0]=16;
	arr[1]=16;
	arr[2]=1;
	arr[3]=(int)((double)NumberOfFrames/FRAMESTEP);
	return 1;
}


public short[] beginframe;
public short[] endframe;
public int setupbaselinesubtract(int start, int end){
//grab frame at beginning and end
  beginframe=new short[256];
  endframe=new short[256];
  shortbuffer.position(start*256);
  shortbuffer.get(beginframe,0,256);
  shortbuffer.position(end*256);
  shortbuffer.get(endframe,0,256);
  baselinesubtract=true;
  firstindex=start;
  lastindex=end;
  return 1;
}


/** routines that can be called directly from matlab **/
public int[] xybfarr;
public int[] matlabReturnXYBandsFrames() {
    if (xybfarr == null) {
	xybfarr = new int[4];
    }
    ReturnXYBandsFrames(xybfarr, 0);
    return xybfarr;
}


/** routines that can be called directly from matlab **/
public int[] framearray;
public int[] matlabUpdateImageArray() {
    if (xybfarr == null) {
	matlabReturnXYBandsFrames();;
    }

    if (framearray == null) {
	framearray = new int[256];
    }

    UpdateImageArray(framearray, xybfarr[0], xybfarr[1], 0);
    return framearray;
}

public int UpdateImageArray(int[] arr, int xdim, int ydim, int instance){

    int val=0;
    try{
    for (int j=0;j<arr.length;j++)arr[j]=0;

    shortbuffer.get(singleframe,0, singleframe.length);
    if (baselinesubtract){
	  for (int j=0;j<arr.length;j++)arr[j]=(int)(singleframe[j]-(beginframe[j]+ (FrameNumber-firstindex)*(((double)(endframe[j]-beginframe[j]))/(lastindex-firstindex))));

	}else
	  for (int j=0;j<arr.length;j++)arr[j]=(int)singleframe[j];

  	FrameNumber++;

   }  catch(Exception e){
	System.out.println("error reading into filebuffer in Hamamatsu Decoder ");
	e.printStackTrace();
}
return 1;
}

/** These get methods aren't used except for testing purposes **/
int getByte(RandomAccessFile f) throws IOException {
		int b = f.read();
		if (b ==-1) throw new IOException("unexpected EOF");
		return b;
	}

/** These get methods aren't used except for testing purposes **/
int getShort(RandomAccessFile f) throws IOException {
	int b0 = getByte(f);
	int b1 = getByte(f);
	return ((b1 << 8) + b0);
}

/** These get methods aren't used except for testing purposes **/
int getHamamatsuShort(RandomAccessFile f) throws IOException{
	int tmp=getShort(f);
	if (tmp>=32768) tmp-=65536;
	return tmp;


}

/** These get methods aren't used except for testing purposes **/
float getFloat(RandomAccessFile f) throws IOException {
	 int accum=0;
	 for (int shiftBy=0;shiftBy<32;shiftBy+=8){
		 accum|=(getByte(f) & 0xff)<<shiftBy;
	 }
	 return Float.intBitsToFloat(accum);


}


/** These get methods aren't used except for testing purposes **/
int getInt(RandomAccessFile f) throws IOException{

	 int accum=0;
		 for (int shiftBy=0;shiftBy<32;shiftBy+=8){
			 accum|=(getByte(f) & 0xff)<<shiftBy;
		 }
		 return accum;

}


/** These get methods aren't used except for testing purposes **/
short getShortInt(RandomAccessFile f) throws IOException{

	 short accum=0;
		 for (int shiftBy=0;shiftBy<16;shiftBy+=8){
			 accum|=(getByte(f) & 0xff)<<shiftBy;
		 }
		 return accum;

}

/** These get methods aren't used except for testing purposes **/
char getChar(RandomAccessFile f) throws IOException{
	int low=getByte(f)&0xff;
	int high=getByte(f);
	return (char)(high<<8|low);

}


public int SumROIs(int[][] rois, String outfile, int startframe, int endframe,int instance){
     System.out.println("in adapter SumROIs start="+startframe+" end="+endframe);
     //determine x,y dimensions
     int[] dims=new int[4];
     ReturnXYBandsFrames(dims, 0);
     int framewidth=dims[0];
     int frameheight=dims[1];
     int numframes=dims[3];
     System.out.println("debug width ="+framewidth+" height ="+frameheight+" num frames ="+numframes);
     //bounds check
     if (endframe<0) endframe=numframes; /*a shortcut for scanning whole record*/
	 if (endframe>numframes) endframe=numframes;
	 if (startframe<0) startframe=0;
     if (startframe>numframes) startframe=0;

     //create an array to hold the sums from the rois.
     int[][] sum=new int[endframe-startframe][rois.length];

     //create an array to read the data into.
     int[] tmparray=new int[framewidth*frameheight];

     //this special case checks for unit rois and treats them differently
     boolean unitroi=true;
     for (int i=0;i<rois.length;i++){
		 if (rois[i].length!=1) unitroi=false;
		 }



	 if (unitroi){
		 System.out.println("USING UNIT ROI Hamamatsu METHOD");
		 for (int k=startframe;k<endframe;k++){
		 for (int i=0;i<rois.length;i++){
		 for (int f=0;f<FRAMESTEP;f++){
		  sum[k-startframe][i]=(int)buffer.get(rois[i][0]+(int)((k*FRAMESTEP+f)*sNumChns));//getHamamatsuShort(Picfile);
	    }
		}
   		}
	 }else{


     for (int k=startframe;k<endframe;k++){//for each frame
     JumpToFrame(k,0); //goto frame
     UpdateImageArray(tmparray,framewidth,frameheight,0); //load image
     for (int i=0;i<rois.length;i++){ //for each roi

       for (int j=0;j<rois[i].length;j++){ //go to each element in the roi
         sum[k-startframe][i]+=tmparray[rois[i][j]];

       }//j
       //System.out.println("sim for frame "+k+" = "+sum[k-startframe][i]);
      }//i
     }//k
	}//else
    //generate output
    try{
    PrintWriter file=new PrintWriter(new FileWriter(outfile),true);
    System.out.println("printing roi file end "+endframe+" start "+startframe);
    for (int j=startframe;j<endframe;j++){
	 // System.out.println("debug frame="+j);
	 file.print((j)+" ");
	// System.out.print("debug "+(j)+" ");
	 for (int i=0;i<rois.length;i++) {
		 //System.out.println("debug rois "+i+" out of "+rois.length);
		 file.print(sum[(j-startframe)][i]+" ");
		 //System.out.print("debug "+sum[(j-startframe)][i]+" ");
		 }
	 file.print("\n");
	 //System.out.print("\n");
	 //System.out.println("debug frame="+j);
	 }
    file.close();
    System.out.println("debug closed file");
    }catch(IOException e){System.out.println("error opening file for rois...");}
     catch(Exception e){System.out.println("Some other error in sumROis");e.printStackTrace();}
    return 1;
   }



public String export_valseparator=" ";
public String export_lineseparator="\n";
public String export_frameseparator="\n\n";

/** Utility method for export frame data to text files. Set global variables
export_valseparator, export_lineseparator, export_frameseparator to control output.
(export_lineseparator and export_frameseparator can be set to null to generate one
long line of data). Set endfame to -1 to have endframe=maxnumber frames.
*/

public void ExportFrames(String outfile,int startframe, int endframe){
	int[] info=new int[4];
	ReturnXYBandsFrames(info, 0);
	int xdim=info[0];
	int ydim=info[1];
	int[] frame=new int[xdim*ydim];
	if (endframe == -1) endframe=info[3];
	if (startframe < 0) startframe=0;
	if (endframe <= startframe) endframe=startframe+1;

	int ystart=0;
	int yend=16;

	try{
	    PrintWriter file=new PrintWriter(new FileWriter(outfile),true);
	    for (int j=startframe;j<endframe;j++){
		 UpdateImageArray(frame,xdim,ydim,0);
		 framescan:
		 for (int y=ystart;y<yend;y++){

		  for (int x=0;x<xdim;x++) {
			 file.print(frame[y*xdim+x]+export_valseparator);
			 }
		  if(export_lineseparator!=null) file.print(export_lineseparator);
		 }
		 if (export_frameseparator!=null) file.print(export_frameseparator);
		 }
	    file.close();

	    }catch(IOException e){System.out.println("error opening file for rois...");}
}

/** Short form for generating a comma seperated list of numbers**/
public void ExportFramesCSL(String outfile, int startframe, int endframe){
	export_valseparator=",";
	export_lineseparator=null;
	export_frameseparator=null;
	ExportFrames(outfile,startframe,endframe);

}










/*
public static void main(String[] arg){
System.out.println("hamamatsu reader test");
HamamatsuDecoder pd=new HamamatsuDecoder();
pd.OpenImageFile(arg[0]);
pd.test();



}
*/

public static void main(String[] args){
System.out.println("ARGUS  reader");
if (args.length!=4){
	System.out.println(
		 "export .PDA files to ascii.  Usage:\n"+
		 "java HamamatsuDecoder infile, outfile, startframenumber, endframenumber \n"+
		 "where infile is a raw data file from Hamamatsu's Argus optical mapping system\n"+
		 "      start and endframenumber give the range of frames to export to ascii\n"+
		 "      (*** where endframenumber==-1 is short for max number of frames)\n"+
		 "      \n");
  return;
}
HamamatsuDecoder pd=new HamamatsuDecoder();

int success=pd.OpenImageFile(args[0]);
if (success==0) {
	System.out.println("Error reading input file");
	return;
}
pd.ExportFrames(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
}


}


