int i = 0;
void setup() {
  //initialize serial communications at a 9600 baud rate
  Serial.begin(9600, SERIAL_8E2);
}

void loop()
{
  //send 'Hello, world!' over the serial port
  i = (i + 1) % 1024;
  Serial.println(i);
  //wait 100 milliseconds so we don't drive ourselves crazy
  delay(10);
}
