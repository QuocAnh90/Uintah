
/*
 *  SoundPort.cc: Handle to the Sound Data type
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <SoundPort.h>
#include <Connection.h>
#include <NotFinished.h>
#include <Port.h>
#include <Classlib/Assert.h>
#include <Classlib/String.h>
#include <iostream.h>

static clString sound_type("Sound");
static clString sound_color("aquamarine4");

SoundIPort::SoundIPort(Module* module, const clString& portname, int protocol)
: IPort(module, sound_type, portname, sound_color, protocol),
  mailbox(10), state(Begin), sample_buf(0)
{
}

SoundIPort::~SoundIPort()
{
    if(sample_buf)
	delete[] sample_buf;
}

void SoundOPort::reset()
{
    state=Begin;
    total_samples=0;
    rate=-1;
}

SoundOPort::SoundOPort(Module* module, const clString& portname, int protocol)
: OPort(module, sound_type, portname, sound_color, protocol),
  sbuf(0), in(0), ptr(0)
{
}

SoundOPort::~SoundOPort()
{
    if(sbuf)
	delete[] sbuf;
}

void SoundIPort::reset()
{
    state=Begin;
}

void SoundIPort::finish()
{
    if(state != Done){
	cerr << "Not all of sound was read...\n";
    }
    if(sample_buf){
	delete[] sample_buf;
	sample_buf=0;
    }
}

int SoundIPort::nsamples()
{
    if(using_protocol() == SoundIPort::Stream){
	cerr << "SoundIPort error: nsamples requested when using Stream protocol\n";
	return 0;
    }
    while(state == Begin)
	do_read();
    return total_samples;
}

double SoundIPort::sample_rate()
{
    while(state == Begin)
	do_read();

    return rate;
}

double SoundIPort::next_sample()
{
    if(state != HaveSamples)
	do_read();
    double s=sample_buf[bufp++];
    if(bufp>=sbufsize){
	state=NeedSamples;
	if(state != HaveSamples)
	    do_read();
	bufp=0;
    }
    return s;
}

int SoundIPort::end_of_stream()
{
    return state==Done;
}

void SoundIPort::do_read()
{
    turn_on();
    SoundComm* comm;
    comm=mailbox.receive();
    switch(comm->action){
    case SoundComm::Parameters:
	ASSERT(state==Begin);
	total_samples=comm->nsamples;
	rate=comm->sample_rate;
	recvd_samples=0;
	state=NeedSamples;
	break;
    case SoundComm::SoundData:
	ASSERT(state==NeedSamples);
	if(sample_buf)
	    delete[] sample_buf;
	sample_buf=comm->samples;
	sbufsize=comm->sbufsize;
	recvd_samples+=sbufsize;
	state=HaveSamples;
	break;
    case SoundComm::EndOfStream:
	ASSERT(state==NeedSamples || state==Begin);
	state=Done;
	total_samples=recvd_samples;
	break;
    }
    delete comm;
    turn_off();
}

void SoundOPort::finish()
{
    // Flush the stream and send an end of stream marker...
    turn_on();
    if(!in){
	Connection* connection=connections[0];
	in=(SoundIPort*)connection->iport;
    }
    if(ptr != 0){
	SoundComm* comm=new SoundComm;
	comm->action=SoundComm::SoundData;
	comm->sbufsize=ptr;
	comm->samples=sbuf;
	in->mailbox.send(comm);
	sbuf=0;
	ptr=0;
    }
    SoundComm* comm=new SoundComm;
    comm->action=SoundComm::EndOfStream;
    in->mailbox.send(comm);
    state=End;
    turn_off();
}

void SoundOPort::set_nsamples(int s)
{
    ASSERT(state == Begin);
    total_samples=s;
}

void SoundOPort::set_sample_rate(double r)
{
    ASSERT(state == Begin);
    rate=r;
}

void SoundOPort::put_sample(double s)
{
    ASSERT(state != End);
    if(state == Begin){
	// Send the Parameters message...
	turn_on();
	SoundComm* comm=new SoundComm;
	comm->action=SoundComm::Parameters;
	comm->sample_rate=rate;
	comm->nsamples=total_samples;
	if(!in){
	    Connection* connection=connections[0];
	    in=(SoundIPort*)connection->iport;
	}
	in->mailbox.send(comm);
	state=Transmitting;
	sbufsize=(int)(rate/20);
	ptr=0;
	turn_off();
    }
    if(!sbuf){
	sbuf=new double[sbufsize];
	ptr=0;
    }
    sbuf[ptr++]=s;
    if(ptr >= sbufsize){
	// Send it away...
	turn_on();
	SoundComm* comm=new SoundComm;
	comm->action=SoundComm::SoundData;
	comm->sbufsize=sbufsize;
	comm->samples=sbuf;
	ASSERT(in != 0);
	in->mailbox.send(comm);
	sbuf=0;
	ptr=0;
	turn_off();
    }
}
