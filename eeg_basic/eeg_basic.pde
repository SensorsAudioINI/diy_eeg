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
import ddf.minim.effects.*;
import ddf.minim.*;

// AUDIO DECLARATIONS
// -------------------------
Minim minim;  
AudioInput line_in;
int instBufferSize   =  1024;
int longerBufferSize = 32768;
float sampleRate     = 22050;
NotchFilter notch;
boolean notchEnabled =  true;
float NOTCH_FREQ     =    50;
float NOTCH_BW       =     1;

// DATA AND FFT DECLARATIONS
// -------------------------
//   FFT objects
FFT fftInst;   // Instantaneous, shorter-buffer FFT
FFT fftLonger; // Longer-buffer FFT that supports lower frequencies
//   Key ranges for EEG
int MAX_FREQ = 100; // The highest frequency we will care about in the zoomed plot
//   Start and end frequencies of the bar, e.g., 0.1-3.5, 3.5-7.5, 7.5-15.5, etc.
static float[] highlightStarts = {0.1, 3.5,  7.5, 15.5, 31.5};
static float[] highlightEnds =   {3.5, 7.5, 15.5, 31.5, 100};
//   Data stores used for calculating average power in a band
float balance = 0.95; // Ratio of running average (95% previous, 5% new)
float [] currPowers = new float[highlightStarts.length];
int   [] currCounts = new   int[highlightStarts.length];
float []  avgPowers = new float[highlightStarts.length];

// DISPLAY DECLARATIONS
// -------------------------
// REDO ALL THIS
float averageScaling = 0.0;
float scalingBalance = 0.90;
float height5; // Height of one of the display bands
static int TEXT_HEIGHT = 25; // Size of text in pixels
static int MAX_HEIGHT = 160-int(TEXT_HEIGHT*1.1); // Drawable size of a band
float fullSpecBalance = 0.99; // Balance for autoscaling - 99% previous, 1% current
float fullSpectrumScale = 0.1; // Spectrum display rescaling factor
float EEGSpecBalance = 0.95; // Balance for autoscaling the 0-100 range: 95% previous, 5% current
float EEGSpectrumScale = 0.1; // Spectrum display rescaling factor
float lineScale = 20;
float lineBalance = 0.95;

// Instantiate the special FFT shift-buffer
FFTBuf longFFTBuf = new FFTBuf(longerBufferSize);

// CODE BEGINS HERE
// --------------------------------------------------------------
void setup()
{
  // Set up display (in pixels)
  size(512, 800);
  // Set up how we draw rectangles (by specifying corners)
  rectMode(CORNERS);  
  // Set up how big each section is
  height5 = height/5;

  // Set up audio
  minim = new Minim(this);
  line_in = minim.getLineIn(Minim.STEREO, instBufferSize, sampleRate);
  line_in.addListener( longFFTBuf );
  notch = new NotchFilter(NOTCH_FREQ, NOTCH_BW, line_in.sampleRate());
  line_in.addEffect( notch );
  
  // Set up normal FFT
  fftInst = new FFT( line_in.bufferSize(), line_in.sampleRate() );
  
  // create an FFT object for calculating slower signals, same sample rate but with more data
  fftLonger = new FFT( longerBufferSize, line_in.sampleRate() );
}

// Render
void draw(){
  
  // Black background
  background(0);
  // Text size is 18pt
  textSize( 18 );
 
  // Display status - the text at the top
  drawTopText();
  
  // Perform a forward FFT on the samples in line_in's mix buffer;
  //   of course, mix vs left vs right doesn't matter if MONO, which it is.
  fftInst.forward( line_in.mix );
  // Also do an FFT of the shift buffer
  fftLonger.forward( longFFTBuf.getSamples() );
  
  // Calculate autoscaling
  calcSpecScales();
  
  // The argument is which row it should occupy
  drawLinePlot(2);
  drawSpectrum(3);
  drawHighResSpectrum(4);
  calcAndDrawEEGBins(5);
}

void calcSpecScales(){
  // Loop over all the spectrogram values, find the maximum power
  float maxVal = 0.0;
  for(int i = 0; i < fftInst.specSize(); i++){
    maxVal = max(maxVal, fftInst.getBand(i));    
  }
  fullSpectrumScale = fullSpecBalance * fullSpectrumScale + (1-fullSpecBalance) * maxVal;
  
  // Loop over the longer FFT buffer, do the same
  float maxLonger = 0.0;
  for(int i = 0; i < fftLonger.specSize(); i++){
    maxLonger = max(maxLonger, fftLonger.getBand(i));    
  }
  EEGSpectrumScale = EEGSpecBalance * EEGSpectrumScale + (1-EEGSpecBalance) * maxLonger;
}

// DRAWING ROUTINES
// --------------------------------------------------------------

// Draw the text at the top
void drawTopText(){
  fill(255, 128); // Color for text
  // Arguments for text is (text, x start, y start)
  text("EEG Spectrum Analyzer", width/4, 0.5*TEXT_HEIGHT);  
  text("FFT res. [q/w] (Hz): " + String.format("%.2f",line_in.sampleRate()/longerBufferSize,3,1), 5, 2.5*TEXT_HEIGHT);
  text("FFT res. [q/w] (s): " + String.format("%.2f",longerBufferSize/line_in.sampleRate(),3,1), 5, 3.5*TEXT_HEIGHT);
  text("Averaging [a/s] (%): " + String.format("%02.2f",balance*100), 5, 4.5*TEXT_HEIGHT);
  text("Notch Filter [z] : " + (notchEnabled? "Yes" : "No"), 5, 5.5*TEXT_HEIGHT);  

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

// Handle user input
void keyPressed()
{
  if ( key == 'q' ){
    longerBufferSize *= 2;
    // Remove the buffer and reinitialize with longer length
    line_in.removeListener( longFFTBuf );
    longFFTBuf = new FFTBuf(longerBufferSize);
    line_in.addListener( longFFTBuf );
    fftLonger = new FFT( longerBufferSize, line_in.sampleRate() );
  }
  if ( key == 'w' ){
    longerBufferSize = max(longerBufferSize/2, line_in.bufferSize());
    // Remove the buffer and reinitialize with shorter length     
    line_in.removeListener( longFFTBuf );
    longFFTBuf = new FFTBuf(longerBufferSize);
    line_in.addListener( longFFTBuf );
    fftLonger = new FFT( longerBufferSize, line_in.sampleRate() );
  }  
  if ( key == 'z' ){
    // Toggle use of the notch filter
    if(notchEnabled){
      line_in.disableEffect( notch );
      notchEnabled = false;
    } else {
      line_in.enableEffect( notch );
      notchEnabled = true;
    }
  }    
  if ( key == 'a' ){
    balance = max(balance-0.005,0);
  } 
  if ( key == 's' ){
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
  colorMode(HSB, highlightStarts.length*5, 255, 255);
  if(stroke){
    stroke(i, 255, 255);
  }
  else{
    fill(i, 255, 255);
  }
  colorMode(RGB, 255, 255, 255, 255);  
}

//   Draw the plot that shows the line over time
void drawLinePlot(int D_IDX){
  // draw the line plot
  noFill();
  int buffer_step = floor(float(longerBufferSize)/width); // Each pixel corresponds to a spot in the buffer
  float [] localFFTSamples = longFFTBuf.getSamples(); // Fetch the samples from the FFT
  float maxCalc = 0;
  for(int i = 0; i < width-1; i++)
  {      
    float brightnessScale = 0.4+0.6*(float(i)/width); // Choose a fading brightness color
    stroke(255*brightnessScale, 255*brightnessScale, 255);      
    float height1 = height5/2 + constrain(height5/2*localFFTSamples[i*buffer_step]/lineScale, -height5/2, height5/2);
    float height2 = height5/2 + constrain(height5/2*localFFTSamples[(i+1)*buffer_step]/lineScale, -height5/2, height5/2);
    line(i, D_IDX*height5 - height1, i+1, D_IDX*height5 - height2);
    maxCalc = max(maxCalc, abs(localFFTSamples[i*buffer_step]), abs(localFFTSamples[(i+1)*buffer_step]));
  }  
  fill(255, 128);
  text("(Autoscaling) Sample History", 5, (D_IDX-1)*height5 + TEXT_HEIGHT);  
  lineScale = lineBalance * lineScale + (1-lineBalance) * maxCalc;  //Autoscale
}

//   Draw the spectrum
void drawSpectrum(int D_IDX){
  noFill();
  float centerFrequency = 0;
  float centerPower = 0;
  for(int i = 0; i < fftInst.specSize(); i++)
  {
    if(getBand(i, fftInst) != -1) {
      setBandColor(getBand(i, fftInst), true);
    }
    else{
      stroke(0, 62, 255);
    }      

    // if the mouse is over the spectrum value we're about to draw
    // set the stroke color to red
    if ( i == mouseX )
    {
      centerFrequency = fftInst.indexToFreq(i);
      centerPower = fftInst.getBand(i);
      stroke(255, 255, 255);
    }
    line(i, D_IDX*height5, i, D_IDX*height5 - min(MAX_HEIGHT, MAX_HEIGHT*fftInst.getBand(i)/fullSpectrumScale));
  }    
  fill(255, 128);
  text("(Autoscaling) Linear Freq: " + String.format("%.1f",centerFrequency) + " Power " + String.format("%.1f",centerPower), 5, (D_IDX-1)*height5 + TEXT_HEIGHT);
}

void drawHighResSpectrum(int D_IDX){
  noStroke();  
  float centerFrequency = 0;
  float centerPower = 0;
  int max_index = fftLonger.freqToIndex(MAX_FREQ);
  int w = ceil( float(width) / (max_index+1) );
  for(int i = 0; i < max_index; i++)
  {
    if(getBand(i, fftLonger) != -1) {
      setBandColor(getBand(i, fftLonger), false);
    }
    else{
      fill(0, 62, 255);        
    }
    // if the mouse is over the spectrum value we're about to draw
    // set the stroke color to red
    if ( mouseX >= i*w && mouseX < i*w + w)
    {
      // Insert label text
      centerFrequency = fftLonger.indexToFreq(i);  
      centerPower = fftLonger.getBand(i);
      // Set up fill for rectangle
      fill(255, 255, 255);        
    }
    rect(i*w, D_IDX*height5, i*w+w, D_IDX*height5 - min(MAX_HEIGHT, MAX_HEIGHT*fftLonger.getBand(i)/EEGSpectrumScale));
  }
  fill(255, 128);    
  text("(Autoscaling) 0-100 Hz Freq: " + String.format("%.1f",centerFrequency) + " Power " + String.format("%.1f",centerPower), 5, (D_IDX-1)*height5 + TEXT_HEIGHT);
}

// Bin the EEG energy and draw it
void calcAndDrawEEGBins(int D_IDX){
  // Reset the power bins
  for(int h_i=0; h_i < currPowers.length; h_i++){
    currPowers[h_i] = 0;
    currCounts[h_i] = 0;
  }
  
  // Assign the FFT power to the appropriate ranges
  int max_index = fftLonger.freqToIndex(highlightEnds[highlightEnds.length-1]);
  for(int i = 0; i <= max_index; i++){
    for(int h_i=0; h_i<highlightStarts.length; h_i++){
      if(fftLonger.freqToIndex(highlightStarts[h_i]) <= i && fftLonger.freqToIndex(highlightEnds[h_i]) > i){
        currPowers[h_i] += fftLonger.getBand(i);
        currCounts[h_i] += 1;
      }
    }   
  }
  
  // Work in the new power averages
  float maxPower = 0;
  for(int i = 0; i < avgPowers.length; i++){     
    if(Float.isNaN(avgPowers[i])){ // Reset if errors
      avgPowers[i] = currPowers[i]/currCounts[i];
    } else {
      avgPowers[i] = avgPowers[i]*balance + (1.0-balance)*currPowers[i]/currCounts[i];
    }    
    maxPower = max(maxPower,avgPowers[i]);    
  }  
  // Adjust the scaling of the bars
  averageScaling = max(maxPower*1, averageScaling);
  averageScaling = averageScaling*scalingBalance + (1-scalingBalance)*maxPower;  
  
  // Redraw the bands
  int w = ceil( float(width) / (highlightStarts.length+1) );  
  int selected = 0;
  for(int h_i=0; h_i < highlightStarts.length; h_i++){
    // if the mouse is over the spectrum value we're about to draw
    // set the stroke color to red
    setBandColor(h_i, false); 
    if ( mouseX >= h_i*w && mouseX < h_i*w + w)
    {
      selected = h_i; // Store which band we've selected        
      fill(255, 255, 255); // Change the highlight fill
    }     
    rect(h_i*w, D_IDX*height5, h_i*w+w, D_IDX*height5 - min(MAX_HEIGHT, avgPowers[h_i]/averageScaling * MAX_HEIGHT));
  }
  fill(255, 128);    
  text("(Autoscaling) Spectrum Selected Frequency: " + highlightStarts[selected] + " - " + highlightEnds[selected], 
    5, (D_IDX-1)*height5 + TEXT_HEIGHT);      
}