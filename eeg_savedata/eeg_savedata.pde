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

// File save
float curr_file_idx = 0;
boolean fileSaving = false;
PrintWriter output;
float decimate = 4;

// AUDIO DECLARATIONS
// -------------------------
Minim minim;  
AudioInput fast_in;
int fastBufferSize = 1024;
int slowBufferSize = 32768;
float fastSampleRate = 22050;
NotchFilter notch;
boolean notchEnabled = true;
float NOTCH_FREQ = 50;
float NOTCH_BW = 1;

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
float balance = 0.95;
float [] currPowers = new float[highlightStarts.length];
int []   currCounts = new int[highlightStarts.length];
float []  avgPowers = new float[highlightStarts.length];

// DISPLAY DECLARATIONS
// -------------------------
float averageScaling = 0.0;
float scalingBalance = 0.90;
float height5; // Height of one of the display bands
static int TEXT_HEIGHT = 25; // Size of text in pixels
static int MAX_HEIGHT = 160-int(TEXT_HEIGHT*1.1); // Size of a band
float specBalance = 0.99; // Balance for autoscaling - 99% previous, 1% current
float spectrumScale = 0.1; // Spectrum display rescaling factor
float highResSpecBalance = 0.95; // Balance for autoscaling - 95% previous, 5% current
float highResSpectrumScale = 0.1; // Spectrum display rescaling factor
float lineScale = 20;
float lineBalance = 0.95;
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
     if(fileSaving && (i%decimate == 0)){
       output.println(newSamples[i]);
     }
    }
    return oldSamples;
  }  
}
FFTBuf slowSamples = new FFTBuf(slowBufferSize);

// CODE BEGINS HERE
// --------------------------------------------------------------
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
  notch = new NotchFilter(NOTCH_FREQ, NOTCH_BW, fastSampleRate);
  fast_in.addEffect( notch );
  
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
  
  // Calculate autoscaling
  calcSpecScales();
  
  // The argument is which row it should occupy
  drawLinePlot(2);
  drawSpectrum(3);
  drawHighResSpectrum(4);
  calcAndDrawEEGBins(5);
}

void calcSpecScales(){
  float maxCalc = 0.0;
  for(int i = 0; i < fftCalc.specSize(); i++){
    maxCalc = max(maxCalc, fftCalc.getBand(i));    
  }
  spectrumScale = specBalance * spectrumScale + (1-specBalance) * maxCalc;
  
  float maxSlow = 0.0;
  for(int i = 0; i < fftSlow.specSize(); i++){
    maxSlow = max(maxSlow, fftSlow.getBand(i));    
  }
  highResSpectrumScale = highResSpecBalance * highResSpectrumScale + (1-highResSpecBalance) * maxSlow;
}

// DRAWING ROUTINES
// --------------------------------------------------------------

// Draw the text at the top
void drawStatus(){
  fill(255, 128); // Color for text
  text("EEG Spectrum Analyzer - Save Data (Press 'M' to Mark)", 20, 0.5*TEXT_HEIGHT);
  text("Saving [k/l] " + (fileSaving ? "(ON): " : "(OFF): ") + String.format("%03.0f",curr_file_idx)+".txt" , 5, 1.5*TEXT_HEIGHT);
  text("FFT res. [q/w] (Hz): " + String.format("%.2f",fast_in.sampleRate()/slowBufferSize,3,1), 5, 2.5*TEXT_HEIGHT);
  text("FFT res. [q/w] (s): " + String.format("%.2f",slowBufferSize/fast_in.sampleRate(),3,1), 5, 3.5*TEXT_HEIGHT);
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
  if ( key == 'z' ){
    if(notchEnabled){
      fast_in.disableEffect( notch );
      notchEnabled = false;
    } else {
      fast_in.enableEffect( notch );
      notchEnabled = true;
    }
  }    
  if ( key == 'a' ){
    balance = max(balance-0.005,0);
  } 
  if ( key == 's' ){
    balance = min(balance+0.005,1);
  }    
  if ( key == 'k' ){
    fileSaving = true;
    output = createWriter("samples_" + String.format("%03.0f",curr_file_idx)+".txt");
    curr_file_idx += 1;
  }      
  if (key == 'l'){
    if(fileSaving){
      output.flush();
      output.close();
      fileSaving = false;
    }
  }
  if(key=='m'){
    output.flush();
    output.println("mark");
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
  int buffer_step = floor(float(slowBufferSize)/width); // Each pixel corresponds to a spot in the buffer
  float [] localSlowSamples = slowSamples.getSamples(); // Fetch the samples from the FFT
  float maxCalc = 0;
  for(int i = 0; i < width-1; i++)
  {      
    float brightnessScale = 0.4+0.6*(float(i)/width); // Choose a fading brightness color
    stroke(255*brightnessScale, 255*brightnessScale, 255);      
    float height1 = height5/2 + constrain(height5/2*localSlowSamples[i*buffer_step]/lineScale, -height5/2, height5/2);
    float height2 = height5/2 + constrain(height5/2*localSlowSamples[(i+1)*buffer_step]/lineScale, -height5/2, height5/2);
    line(i, D_IDX*height5 - height1, i+1, D_IDX*height5 - height2);
    maxCalc = max(maxCalc, abs(localSlowSamples[i*buffer_step]), abs(localSlowSamples[(i+1)*buffer_step]));
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
      centerPower = fftCalc.getBand(i);
      stroke(255, 255, 255);
    }
    line(i, D_IDX*height5, i, D_IDX*height5 - min(MAX_HEIGHT, MAX_HEIGHT*fftCalc.getBand(i)/spectrumScale));
  }    
  fill(255, 128);
  text("(Autoscaling) Linear Freq: " + String.format("%.1f",centerFrequency) + " Power " + String.format("%.1f",centerPower), 5, (D_IDX-1)*height5 + TEXT_HEIGHT);
}

void drawHighResSpectrum(int D_IDX){
  noStroke();  
  float centerFrequency = 0;
  float centerPower = 0;
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
      centerPower = fftSlow.getBand(i);
      // Set up fill for rectangle
      fill(255, 255, 255);        
    }
    rect(i*w, D_IDX*height5, i*w+w, D_IDX*height5 - min(MAX_HEIGHT, MAX_HEIGHT*fftSlow.getBand(i)/highResSpectrumScale));
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
  int max_index = fftSlow.freqToIndex(highlightEnds[highlightEnds.length-1]);
  for(int i = 0; i <= max_index; i++){
    for(int h_i=0; h_i<highlightStarts.length; h_i++){
      if(fftSlow.freqToIndex(highlightStarts[h_i]) <= i && fftSlow.freqToIndex(highlightEnds[h_i]) > i){
        currPowers[h_i] += fftSlow.getBand(i);
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