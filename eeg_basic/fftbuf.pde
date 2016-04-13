// Special threadsafe (for audio) class for shifting in new samples
class FFTBuf implements AudioListener
{
  private float[] left;
  private float[] right;
  private int decimate = 4;
  // We make this public and not thread-safe, but that's okay for now.
  public boolean fileSaving;  
  public PrintWriter output;
  
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