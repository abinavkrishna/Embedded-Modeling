// play.cpp: conceptual audio/video stream decoding example
// Author: RISC project members, CECS
//
// 06/05/18 RD	confirmed clean sc_time constants with GCC 6.3.1
// 04/27/18 RD  clean sc_time constants for newer compilers (e.g. GCC 4.9.2)
// 04/26/18 RD  merged tracing version into main file (define TRACE_FILE)
// 04/24/18 DM  added tracing of stream data (sc_trace is now supported)
// 07/10/17 RD  split the delay into two parts (as suggested by ZC)
// 04/05/17 RD  removed sc_stop from stimulus, let monitor terminate;
//              added COUNT_FRAMES to ensure proper completion
// 03/08/17 RD  avoid duplicated DELAY macros, set BUFSIZE to 10
// 01/13/17 GL  separate version play_fifo.cpp with sc_fifo channels
// 10/04/16 RD  removed obsolete C artifacts from the user-defined queue channel
// 09/30/16 RD  cleanup for release 0.3.0
// 03/07/15 RD  bug fix: use sc_stop() to quit cleanly (exit() is not MT-safe)
// 10/30/15 RD  polished for release 0.2.1
// 09/04/15 GL  use 'long long' (not 'double') to avoid floating-point overflow
//              (FPU exception handling shows high performance variation;
//              moreover, integer arithmetic is more realistic for decoder)
// 08/30/15 GL  initial version based on SpecC example play.sc by RD

#include <stdio.h>
#include "systemc.h"

// frame rates
#define AUDIO_DELAY     26120, SC_US
#define VIDEO_DELAY     33330, SC_US
#define HALF_AUDIO_DELAY 26120/2, SC_US
#define HALF_VIDEO_DELAY 33330/2, SC_US

// simulated time
#define CLIP_LENGTH     10, SC_SEC

// simulated workload unit
#define AUDIOLOAD        7500000
#define VIDEOLOAD       10000000

// illustrate and check frame decoding
//#define PRINT_TIME
//#define PRINT_FRAMES
//#define COUNT_FRAMES
//#define TRACE_FILE

// buffer size in queue channels
#define BUFSIZE		10ul

SC_MODULE(Stimulus)
{
  sc_fifo_out<long long> VideoStream;
  sc_fifo_out<long long> AudioStreamL;
  sc_fifo_out<long long> AudioStreamR;

  void main(void)
  {
    sc_time     t1(0, SC_US);
    sc_time     t2(0, SC_US);
    sc_time	ClipLength(CLIP_LENGTH);
    sc_time	VideoDelay(VIDEO_DELAY);
    sc_time	AudioDelay(AUDIO_DELAY);
    long long   video = 3;
    long long   audioLeft = 27;
    long long   audioRight = 15;

    while ((t1 < ClipLength) || (t2 < ClipLength))
    {
      if (t1 < t2)
      {
        VideoStream.write(video);
        t1 += VideoDelay;
      }
      else
      {
        AudioStreamL.write(audioLeft);
        AudioStreamR.write(audioRight);
        t2 += AudioDelay;
      }
    }
  }

  SC_CTOR(Stimulus)
  {
    SC_THREAD(main);
  }
};

SC_MODULE(VideoCodec)
{
  sc_fifo_in<long long>  p1;
  sc_fifo_out<long long> p2;

  long long d;

  long long decode(long long x)
  {
    int i;

    for (i = 0; i < VIDEOLOAD; i++)
    {
      d *= x;
    }
    return d;
  }

  void main(void)
  {
    long long inFrame, outFrame;

    d = 1;

    while(1)
    {
      inFrame = p1.read();
      wait(HALF_VIDEO_DELAY);
      outFrame = decode(inFrame);
      wait(HALF_VIDEO_DELAY);
      p2.write(outFrame);
    }
  }

  SC_CTOR(VideoCodec)
  {
    SC_THREAD(main);
  }
};

SC_MODULE(AudioCodec)
{
  sc_fifo_in<long long>  p1;
  sc_fifo_out<long long> p2;

  long long d;

  long long decode(long long x)
  {
    int i;

    for (i = 0; i < AUDIOLOAD; i++)
    {
      d *= x;
    }
    return d;
  }

  void main(void)
  {
    long long inFrame, outFrame;

    d = 1;

    while(1)
    {
      inFrame = p1.read();
      wait(HALF_AUDIO_DELAY);
      outFrame = decode(inFrame);
      wait(HALF_AUDIO_DELAY);
      p2.write(outFrame);
    }
  }

  SC_CTOR(AudioCodec)
  {
    SC_THREAD(main);
  }
};

SC_MODULE(DUT)
{
  sc_fifo_in<long long>  VideoIn;
  sc_fifo_in<long long>  AudioInL;
  sc_fifo_in<long long>  AudioInR;
  sc_fifo_out<long long> VideoOut;
  sc_fifo_out<long long> AudioOutLeft;
  sc_fifo_out<long long> AudioOutRight;

  VideoCodec video;
  AudioCodec left, right;

  void before_end_of_elaboration()
  {
    video.p1.bind(VideoIn);
    video.p2.bind(VideoOut);

    left.p1.bind(AudioInL);
    left.p2.bind(AudioOutLeft);

    right.p1.bind(AudioInR);
    right.p2.bind(AudioOutRight);
  }

  SC_CTOR(DUT): video("video"), left("left"), right("right")
  { }
};

SC_MODULE(Display)
{
  sc_fifo_in<long long> VideoIn;
#ifdef COUNT_FRAMES
  unsigned FramesDisplayed;
#endif

  void display(long long /* Frame */)
  {
#ifdef PRINT_FRAMES
    putchar('V');       // display the frame
    fflush(stdout);
#endif
#ifdef COUNT_FRAMES
    FramesDisplayed++;
#endif
  }

  void main(void)
  {
    long long Frame;

    while(1)
    {
      Frame = VideoIn.read();
      display(Frame);
    }
  }

  SC_CTOR(Display)
  {
#ifdef COUNT_FRAMES
    FramesDisplayed = 0;
#endif
    SC_THREAD(main);
  }
};

SC_MODULE(Speaker)
{
  sc_fifo_in<long long> AudioIn;
#ifdef COUNT_FRAMES
  unsigned FramesPlayed;
#endif

  char Channel;

  void sound(long long /* Frame */)
  {
#ifdef PRINT_FRAMES
    putchar(Channel); // output the sound
    fflush(stdout);
#endif
#ifdef COUNT_FRAMES
    FramesPlayed++;
#endif
  }

  void main(void)
  {
    long long Frame;

    while(1)
    {
      Frame = AudioIn.read();
      sound(Frame);
    }
  }

  Speaker(sc_module_name name, char chn): sc_module(name), Channel(chn)
  {
#ifdef COUNT_FRAMES
    FramesPlayed = 0;
#endif
    SC_THREAD(main);
  }
  SC_HAS_PROCESS(Speaker);
};

SC_MODULE(Monitor)
{
  sc_fifo_in<long long> Video;
  sc_fifo_in<long long> AudioLeft;
  sc_fifo_in<long long> AudioRight;

  Display screen;
  Speaker speakerL, speakerR;

  void before_end_of_elaboration()
  {
    screen.VideoIn.bind(Video);
    speakerL.AudioIn.bind(AudioLeft);
    speakerR.AudioIn.bind(AudioRight);
  }

  SC_CTOR(Monitor):
        screen("screen"),
        speakerL("speakerL", 'L'),
        speakerR("speakerR", 'R')
  { }
};

SC_MODULE(Top)
{
  sc_fifo<long long> VideoStream;
  sc_fifo<long long> AudioStreamL;
  sc_fifo<long long> AudioStreamR;
  sc_fifo<long long> Video;
  sc_fifo<long long> AudioLeft;
  sc_fifo<long long> AudioRight;

  Stimulus       stimulus;
  DUT            dut;
  Monitor        monitor;

  void before_end_of_elaboration()
  {
    stimulus.VideoStream.bind(VideoStream);
    stimulus.AudioStreamL.bind(AudioStreamL);
    stimulus.AudioStreamR.bind(AudioStreamR);

    dut.VideoIn.bind(VideoStream);
    dut.AudioInL.bind(AudioStreamL);
    dut.AudioInR.bind(AudioStreamR);
    dut.VideoOut.bind(Video);
    dut.AudioOutLeft.bind(AudioLeft);
    dut.AudioOutRight.bind(AudioRight);

    monitor.Video.bind(Video);
    monitor.AudioLeft.bind(AudioLeft);
    monitor.AudioRight.bind(AudioRight);
  }

  SC_CTOR(Top)
  : VideoStream("VideoStream", BUFSIZE)
  , AudioStreamL("AudioStreamL", BUFSIZE)
  , AudioStreamR("AudioStreamR", BUFSIZE)
  , Video("Video", BUFSIZE)
  , AudioLeft("AudioLeft", BUFSIZE)
  , AudioRight("AudioRight", BUFSIZE)
  , stimulus("stimulus")
  , dut("dut")
  , monitor("monitor")
  { }
};

int sc_main(int argc, char* argv[])
{
  Top top("top");
  sc_time	ClipLength(CLIP_LENGTH);
  sc_time	VideoDelay(VIDEO_DELAY);
  sc_time	AudioDelay(AUDIO_DELAY);
#ifdef TRACE_FILE
  sc_trace_file *f = sc_create_vcd_trace_file(argv[0]);
  f->set_time_unit(sc_time(10, SC_US).to_seconds(), SC_SEC);
  sc_trace(f, top.dut.video.d, "top.dut.video.d");
  sc_trace(f, top.dut.left.d, "top.dut.left.d");
  sc_trace(f, top.dut.right.d, "top.dut.right.d");
#  ifdef COUNT_FRAMES
  sc_trace(f, top.monitor.screen.FramesDisplayed, "top.monitor.screen.FramesDisplayed");
  sc_trace(f, top.monitor.speakerL.FramesPlayed, "top.monitor.speakerL.FramesPlayed");
  sc_trace(f, top.monitor.speakerR.FramesPlayed, "top.monitor.speakerR.FramesPlayed");
#  endif
#endif

#ifdef PRINT_TIME
  printf("%s: ", sc_time_stamp().to_string().c_str());
#endif
  printf("Playing...\n");
  sc_start();
#ifdef PRINT_FRAMES
  printf("\n");
#endif
#ifdef COUNT_FRAMES
  printf("%u video frames displayed\n", top.monitor.screen.FramesDisplayed);
  printf("%u left audio frames played\n", top.monitor.speakerL.FramesPlayed);
  printf("%u right audio frames played\n", top.monitor.speakerR.FramesPlayed);
  assert(top.monitor.screen.FramesDisplayed ==
		ClipLength.value()/VideoDelay.value()+1);
  assert(top.monitor.speakerL.FramesPlayed ==
		ClipLength.value()/AudioDelay.value()+1);
  assert(top.monitor.speakerR.FramesPlayed ==
		ClipLength.value()/AudioDelay.value()+1);
#endif
#ifdef PRINT_TIME
  printf("%s: ", sc_time_stamp().to_string().c_str());
  assert(sc_time_stamp() >= ClipLength);
#endif
#ifdef TRACE_FILE
  sc_close_vcd_trace_file(f);
#endif
  printf("Done!\n");
  return 0;
}

// EOF
