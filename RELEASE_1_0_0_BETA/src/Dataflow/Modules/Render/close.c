/*
  The contents of this file are subject to the University of Utah Public
  License (the "License"); you may not use this file except in compliance
  with the License.
  
  Software distributed under the License is distributed on an "AS IS"
  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
  License for the specific language governing rights and limitations under
  the License.
  
  The Original Source Code is SCIRun, released March 12, 2001.
  
  The Original Source Code was developed by the University of Utah.
  Portions created by UNIVERSITY are Copyright (C) 2001, 1994 
  University of Utah. All Rights Reserved.
*/

/*
 *	iclose and iflush -
 *
 *				Paul Haeberli - 1984
 *
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	"imagelib.h"

int iclose(IMAGE *image)
{
    long tablesize, ret;

    iflush(image);
    img_optseek(image, 0);
    if (image->flags&_IOWRT) {
	if(image->dorev)
	    cvtimage(image);
	if (img_write(image,(char *)image,sizeof(IMAGE)) != sizeof(IMAGE)) {
	    i_errhdlr("iclose: error on write of image header\n");
	    return EOF;
	}
	if(image->dorev)
	    cvtimage(image);
	if(ISRLE(image->type)) {
	    img_optseek(image, 512L);
	    tablesize = image->ysize*image->zsize*sizeof(long);
	    if(image->dorev)
		cvtlongs(image->rowstart,tablesize);
	    if (img_write(image,(char *)(image->rowstart),tablesize) != tablesize) {
		i_errhdlr("iclose: error on write of rowstart\n");
		return EOF;
	    }
	    if(image->dorev)
		cvtlongs(image->rowsize,tablesize);
	    if (img_write(image,(char *)(image->rowsize),tablesize) != tablesize) {
		i_errhdlr("iclose: error on write of rowsize\n");
		return EOF;
	    }
	}
    }
    if(image->base) {
	free(image->base);
	image->base = 0;
    }
    if(image->tmpbuf) {
	free(image->tmpbuf);
	image->tmpbuf = 0;
    }
    if(ISRLE(image->type)) {
	free(image->rowstart);
	image->rowstart = 0;
	free(image->rowsize);
	image->rowsize = 0;
    }
    ret = close(image->file);
    if(ret != 0) 
	i_errhdlr("iclose: error on close of file\n");
    free(image);
    return ret;
}

int iflush(IMAGE *image)
{
    unsigned short *base;

    if ( (image->flags&_IOWRT)
     && (base=image->base)!=NULL && (image->ptr-base)>0) {
	    if (putrow(image, base, image->y,image->z)!=image->xsize) {
		    image->flags |= _IOERR;
		    return(EOF);
	    }
    }
    return(0);
}
