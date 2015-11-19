/**
  * 
  * This code is used to calculate the FFTs of an input signal, with a special focus on frequency ranges that
  * might be relevant for EEGs.  We collect samples frequently, and shift them into an audio buffer that calculates
  * a longer timelength (and thus higher resolution) FFT signal.  Resolution can be changed on the fly.  The only
  * tricky bits is that we register our buffer as a Listener to audio events so that each buffer capture shifts in
  * the data in a thread-safe asynchronous way.  Rest of the code is basically for display.
  *
  * Email daniel.l.neil@gmail.com for comments / issues / questions.
  *
  */

// Add imports for audio manipulation
import ddf.minim.analysis.*;
import ddf.minim.*;

// AUDIO DECLARATIONS
// -------------------------
Minim minim;  
AudioInput fast_in;
int fastBufferSize = 1024;
int slowBufferSize = 8192;
float fastSampleRate = 22050;

// DATA AND FFT DECLARATIONS
// -------------------------
//   FFT objects
FFT fftCalc;
FFT fftSlow;
//   Key ranges for EEG
int MAX_FREQ = 100; // The highest frequency we will care about in the zoomed plot
static float[] highlightStarts = {0.1, 3.5,  7.5, 15.5, 31.5};
static float[] highlightEnds =   {3.5, 7.5, 15.5, 31.5, 100};
//   Data stores used for calculating average power in a band
float balance = 0.75;
float [] currPowers = new float[highlightStarts.length];
int []   currCounts = new int[highlightStarts.length];
float []  avgPowers = new float[highlightStarts.length];

// DISPLAY DECLARATIONS
// -------------------------
float height5; // Height of one of the display bands
static int TEXT_HEIGHT = 25; // Size of text in pixels
static int MAX_HEIGHT = 160; // Size of a band
static float spectrumScale = 4; // Spectrum display rescaling factor

// Special threadsafe class for shifting in new samples
class FFTBuf implements AudioListener
{
  private float[] left;
  private float[] right;
  
  FFTBuf(int num_samples){
    left = new float[num_samples];
    right = new float[num_samples];
    for(int i=0; i<num_samples;i++){
      left[i] = 0;
      right[i] = 0;
    }    
  }
  public synchronized void samples(float[] samp){
    left = addNewSamples(left, samp);
  }
  
  public synchronized void samples(float[] sampL, float[] sampR){
    left = addNewSamples(left, sampL);
    right = addNewSamples(right, sampR);
  }
  
  public synchronized float[] getSamples(){
    return left;
  }
  
  private float[] addNewSamples(float[] oldSamples, float [] newSamples){  
    // Shift buffer left
    for(int i=0; i < oldSamples.length - newSamples.length; i++){
     oldSamples[i] = oldSamples[i+newSamples.length];   
    }  
    // Add in new samples
    for(int i=0; i < newSamples.length; i++){
     oldSamples[i+(oldSamples.length-newSamples.length)] = newSamples[i];   
    }
    return oldSamples;
  }  
}
FFTBuf slowSamples = new FFTBuf(slowBufferSize);

// CODE BEGINS HERE
// -------------------------
void setup()
{
  // Set up display
  size(512, 800);
  rectMode(CORNERS);  
  height5 = height/5;

  // Set up audio
  minim = new Minim(this);
  fast_in = minim.getLineIn(Minim.STEREO, fastBufferSize, fastSampleRate);
  fast_in.addListener( slowSamples );

  // Set up normal FFT
  fftCalc = new FFT( fast_in.bufferSize(), fast_in.sampleRate() );
  
  // create an FFT object for calculating slower signals, same sample rate but with more data
  fftSlow = new FFT( slowBufferSize, fast_in.sampleRate() );
}

// Render
void draw(){
  
  background(0);
  textSize( 18 );
 
  // Display status
  drawStatus();
  
  // perform a forward FFT on the samples in in's mix buffer
  // note that if in were a MONO file, this would be the same as using in.left or in.right
  fftCalc.forward( fast_in.mix );
  fftSlow.forward( slowSamples.getSamples() );
  
  // The argument is which row it should occupy
  drawLinePlot(2);
  drawSpectrum(3);
  drawHighResSpectrum(4);
  calcAndDrawEEGBins(5);
}


// Draw the text at the top
void drawStatus(){
  fill(255, 128); // Color for text
  text("EEG Spectrum Analyzer", width/4, 0.5*TEXT_HEIGHT);  
  text("FFT res. [q/w] (Hz): " + String.format("%.2f",fast_in.sampleRate()/slowBufferSize,3,1), 5, 2.5*TEXT_HEIGHT);
  text("FFT res. [q/w] (s): " + String.format("%.2f",slowBufferSize/fast_in.sampleRate(),3,1), 5, 3.5*TEXT_HEIGHT);
  text("Max Freq. [a/s] (Hz): " + MAX_FREQ, 5, 4.5*TEXT_HEIGHT);
  text("Averaging [z/x] (%): " + String.format("%02.2f",balance*100), 5, 5.5*TEXT_HEIGHT);

  // Right column
  text("Delta (0.1-3 Hz): "  + String.format("%06.2f",avgPowers[0]), width/2, 1.5*TEXT_HEIGHT);
  text("Theta (4-7 Hz): "    + String.format("%06.2f",avgPowers[1]), width/2, 2.5*TEXT_HEIGHT);
  text("Alpha (8-15 Hz): "   + String.format("%06.2f",avgPowers[2]), width/2, 3.5*TEXT_HEIGHT);
  text("Beta (16-31 Hz): "   + String.format("%06.2f",avgPowers[3]), width/2, 4.5*TEXT_HEIGHT);
  text("Gamma (32-100 Hz): " + String.format("%06.2f",avgPowers[4]), width/2, 5.5*TEXT_HEIGHT);
  
  // Draw border
  stroke(128,128,128);
  rect(0, height5-2, width, height5);
}

// Handle changes
void keyPressed()
{
  if ( key == 'q' ){
    slowBufferSize *= 2;
    fast_in.removeListener( slowSamples );
    slowSamples = new FFTBuf(slowBufferSize);
    fast_in.addListener( slowSamples );
    fftSlow = new FFT( slowBufferSize, fast_in.sampleRate() );
  }
  if ( key == 'w' ){
    slowBufferSize = max(slowBufferSize/2, fast_in.bufferSize());
    fast_in.removeListener( slowSamples );
    slowSamples = new FFTBuf(slowBufferSize);
    fast_in.addListener( slowSamples );
    fftSlow = new FFT( slowBufferSize, fast_in.sampleRate() );
  }  
  if ( key == 'a' ){
    MAX_FREQ -=1;
  }    
  if ( key == 's' ){
    MAX_FREQ +=1;
  }   
  if ( key == 'z' ){
    balance = max(balance-0.005,0);
  } 
  if ( key == 'x' ){
    balance = min(balance+0.005,1);
  }    
}

// DISPLAY ROUTINES
// ----------------------------------------------------------------
//   Find which frequency band an index belongs to, or -1 if none.
int getBand(int index, FFT myFFT){
  for(int i=0; i<highlightStarts.length; i++){
    if(myFFT.freqToIndex(highlightStarts[i]) <= index && myFFT.freqToIndex(highlightEnds[i]) > index){
      return i;
    }
  }
  return -1;
}

//   Set the color for a band
void setBandColor(int i, boolean stroke){
  colorMode(HSB, highlightStarts.length*8);
  if(stroke){
    stroke(i, 100, 100);
  }
  else{
    fill(i, 100, 100);
  }
  colorMode(RGB, 255, 255, 255, 255);  
}

//   Draw the plot that shows the line over time
void drawLinePlot(int D_IDX){
  // draw the line plot
  noFill();
  int buffer_step = floor(float(slowBufferSize)/width); // Each pixel corresponds to a spot in the buffer
  float [] localSlowSamples = slowSamples.getSamples(); // Fetch the samples from the FFT
  for(int i = 0; i < width-1; i++)
  {      
    float brightnessScale = 0.25+0.75*(float(i)/width);
    stroke(0, 62*brightnessScale, 255*brightnessScale);        
    line(i, D_IDX*height5 - height5/2 + height5/2*localSlowSamples[i*buffer_step],
         i+1, D_IDX*height5 - height5/2 + height5/2*localSlowSamples[(i+1)*buffer_step]);           
  }  
}

//   Draw the spectrum
void drawSpectrum(int D_IDX){
  noFill();
  float centerFrequency = 0;
  for(int i = 0; i < fftCalc.specSize(); i++)
  {
    if(getBand(i, fftCalc) != -1) {
      setBandColor(getBand(i, fftCalc), true);
    }
    else{
      stroke(0, 62, 255);
    }      

    // if the mouse is over the spectrum value we're about to draw
    // set the stroke color to red
    if ( i == mouseX )
    {
      centerFrequency = fftCalc.indexToFreq(i);
      stroke(255, 255, 255);
    }
    line(i, D_IDX*height5, i, D_IDX*height5 - min(MAX_HEIGHT, fftCalc.getBand(i)*spectrumScale));
  }    
  fill(255, 128);
  text("Spectrum Selected Frequency: " + centerFrequency, 5, D_IDX*height5 - TEXT_HEIGHT);
}

void drawHighResSpectrum(int D_IDX){
  noStroke();  
  float centerFrequency = 0;
  int max_index = fftSlow.freqToIndex(MAX_FREQ);
  int w = ceil( float(width) / (max_index+1) );
  for(int i = 0; i < max_index; i++)
  {
    if(getBand(i, fftSlow) != -1) {
      setBandColor(getBand(i, fftSlow), false);
    }
    else{
      fill(0, 62, 255);        
    }
    // if the mouse is over the spectrum value we're about to draw
    // set the stroke color to red
    if ( mouseX >= i*w && mouseX < i*w + w)
    {
      // Insert label text
      centerFrequency = fftSlow.indexToFreq(i);        
      // Set up fill for rectangle
      fill(255, 255, 255);        
    }
    rect(i*w, D_IDX*height5, i*w+w, D_IDX*height5 - min(MAX_HEIGHT, fftSlow.getBand(i)*spectrumScale));
  }
  fill(255, 128);    
  text("Spectrum Selected Frequency: " + centerFrequency, 5, 4*height5 - TEXT_HEIGHT);    
}

// Bin the EEG energy and draw it
void calcAndDrawEEGBins(int D_IDX){
  // Reset the power bins
  for(int h_i=0; h_i < currPowers.length; h_i++){
    currPowers[h_i] = 0;
    currCounts[h_i] = 0;
  }
  
  // Assign the FFT power to the appropriate ranges
  int max_index = fftSlow.freqToIndex(highlightEnds[highlightEnds.length-1]);
  for(int i = 0; i <= max_index; i++){
    for(int h_i=0; h_i<highlightStarts.length; h_i++){
      if(fftSlow.freqToIndex(highlightStarts[h_i]) <= i && fftSlow.freqToIndex(highlightEnds[h_i]) > i){
        currPowers[h_i] += fftSlow.getBand(i)*spectrumScale;
        currCounts[h_i] += 1;
      }
    }
  }
  
  // Work in the new power averages
  for(int i = 0; i < avgPowers.length; i++){    
    if(Float.isNaN(avgPowers[i])){
      avgPowers[i] = currPowers[i]/currCounts[i];
    } else {
      avgPowers[i] = avgPowers[i]*balance + (1.0-balance)*currPowers[i]/currCounts[i];
    }    
  }  
    
  // Redraw the bands
  int w = ceil( float(width) / (highlightStarts.length+1) );  
  int selected = 0;
  for(int h_i=0; h_i < highlightStarts.length; h_i++){
    avgPowers[h_i] /= currCounts[h_i];
    // if the mouse is over the spectrum value we're about to draw
    // set the stroke color to red
    setBandColor(h_i, false); 
    if ( mouseX >= h_i*w && mouseX < h_i*w + w)
    {
      selected = h_i; // Store which band we've selected        
      fill(255, 255, 255); // Change the highlight fill
    }     
    rect(h_i*w, D_IDX*height5, h_i*w+w, D_IDX*height5 - min(MAX_HEIGHT, avgPowers[h_i]));
  }
  fill(255, 128);    
  text("Spectrum Selected Frequency: " + highlightStarts[selected] + " - " + highlightEnds[selected], 
    5, D_IDX*height5 - TEXT_HEIGHT);      
}