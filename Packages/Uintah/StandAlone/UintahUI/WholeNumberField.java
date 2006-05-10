//**************************************************************************
// Program : WholeNumberField.java
// Purpose : An extension of JTextField to take care of whole numbers
// Author  : Biswajit Banerjee
// Date    : 12/7/1998
// Mods    :
//**************************************************************************
// $Id: WholeNumberField.java,v 1.2 2000/02/03 05:36:59 bbanerje Exp $

//************ IMPORTS **************
import java.awt.*;
import java.text.NumberFormat;
import java.text.ParseException;
import javax.swing.*;
import javax.swing.text.*;

//**************************************************************************
// Class   : WholeNumberField
// Purpose : Creates a text field that validates whole numbers
//**************************************************************************
public class WholeNumberField extends JTextField {

  // Data
  private Toolkit toolkit;
  private NumberFormat integerFormatter;

  // Data that may be needed later

  // Constructor
  public WholeNumberField(int value, int columns) {
    
    // Set the size of the component
    super(columns);
    
    // Get the toolkit
    toolkit = Toolkit.getDefaultToolkit();
    integerFormatter = NumberFormat.getNumberInstance();
    integerFormatter.setParseIntegerOnly(true);
    setValue(value);
  }

  // Get method
  public int getValue() {
    int retVal = 0;
    try {
      retVal = integerFormatter.parse(getText()).intValue();
    } catch (ParseException e) {
    }
    return retVal;
  }

  // set method
  public void setValue(int value) {
    setText(integerFormatter.format(value));
  }

  // Create the related document
  protected Document createDefaultModel() {
    return new WholeNumberDocument();
  }

  // Inner class for whole number document
  protected class WholeNumberDocument extends PlainDocument {

    // The insert string method
    public void insertString(int offs, String src, AttributeSet a)
      throws BadLocationException {

      char[] source = src.toCharArray();
      char[] result = new char[source.length];
      int j = 0;

      for (int i = 0; i < result.length; i++) {
	if (Character.isDigit(source[i])) result[j++] = source[i];
      }
      super.insertString(offs, new String(result,0,j), a);
    }
  }
}
// $Log: WholeNumberField.java,v $
// Revision 1.2  2000/02/03 05:36:59  bbanerje
// Just a few changes in all the java files .. and some changes in
// GenerateParticleFrame and Particle and ParticleList
//
