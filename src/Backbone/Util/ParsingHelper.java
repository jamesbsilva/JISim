package Backbone.Util;
/**
* 
*    @(#) ParsingHelper
*/ 
import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** 
*   A parsing helping class.
*  <br>
* 
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2012-07    
*/
public class ParsingHelper {
    String testStr= "  telte    :temperature:  3845  false terj 9384.9  \\ 938475";
    String testLoc= "/home/j2/JISim/Settings/TestParam.txt";
    String delim=":";
   
    /**
    *   checkParameterExist checks that the parameter exist in the input file.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @return true - if parameter exists
    */
    public boolean checkParameterExist(String filename, String parameter){
        //System.out.println("CHECK FOR : "+parameter);
        if(!parseFileForLine(filename, parameter).equals("PARAMNOTFOUND")){
            return true;
        }else{return false;}
    }
    
    /**
    *   checkParameterExist checks that the parameter exist in the input file.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @return true - if parameter exists
    */
    public boolean checkParameterExist(ArrayList<String> filename, String parameter){
        //System.out.println("CHECK FOR : "+parameter);
        if(!parseStringsForLine(filename, parameter).equals("PARAMNOTFOUND")){
            return true;
        }else{return false;}
    }
    
    /**
    *        parseFileDouble parses a file based on the parameter and returns the double 
    *   found.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return double parsed
    */
    public double parseFileDouble(String filename, String parameter,int locationOffset){
        return parseDouble(parseFileForLine(filename,parameter),parameter,locationOffset);
    }
    
    /**
    *        parseFileDouble parses a file in string list form based on the parameter and returns the double 
    *   found.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return double parsed
    */
    public double parseFileDouble(ArrayList<String> filename, String parameter,int locationOffset){
        return parseDouble(parseStringsForLine(filename,parameter),parameter,locationOffset);
    }
    /**
    *        parseFileInt parses a file based on the parameter and returns the integer 
    *   found.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return integer parsed
    */
    public int parseFileInt(String filename, String parameter,int locationOffset){
        return parseInt(parseFileForLine(filename,parameter),parameter,locationOffset);
    }

    /**
    *        parseFileInt parses a file in string list form  based on the parameter and returns the integer 
    *   found.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return integer parsed
    */
    public int parseFileInt(ArrayList<String> filename, String parameter,int locationOffset){
        return parseInt(parseStringsForLine(filename,parameter),parameter,locationOffset);
    }
    
    /**
    *        parseFileBoolean parses a file based on the parameter and returns the boolean 
    *   found.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return boolean parsed
    */
    public boolean parseFileBoolean(String filename, String parameter,int locationOffset){
        return parseBoolean(parseFileForLine(filename,parameter),parameter,locationOffset);
    }

        /**
    *        parseFileBoolean parses a file in string list form based on the parameter and returns the boolean 
    *   found.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return boolean parsed
    */
    public boolean parseFileBoolean(ArrayList<String> filename, String parameter,int locationOffset){
        return parseBoolean(parseStringsForLine(filename,parameter),parameter,locationOffset);
    }

    
    /**
    *        parseFileFirstString parses a file based on the parameter and returns the first string of data 
    *   found.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return string parsed
    */
    public String parseFileFirstString(String filename, String parameter){
        return parse(parseFileForLine(filename,parameter),parameter,0);
    }

        /**
    *        parseFileFirstString parses a file in string list form based on the parameter and returns the first string of data 
    *   found.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return string parsed
    */
    public String parseFileFirstString(ArrayList<String> filename, String parameter){
        return parse(parseStringsForLine(filename,parameter),parameter,0);
    }
    
    /**
    *        parseFileString parses a file based on the parameter and returns the  string of data 
    *   found in queried column.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return string parsed
    */    
    public String parseFileString(String filename, String parameter, int locationOffset){
        return parse(parseFileForLine(filename,parameter),parameter,locationOffset);
    }

        /**
    *        parseFileString parses a file in string list form based on the parameter and returns the  string of data 
    *   found in queried column.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return string parsed
    */    
    public String parseFileString(ArrayList<String> filename, String parameter, int locationOffset){
        return parse(parseStringsForLine(filename,parameter),parameter,locationOffset);
    }
    
    /**
    *        parseFileForLine parses a file based on the parameter and returns the line 
    *   containing the parameter.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @return string parsed
    */
    public String parseFileForLine(String file,String parameter){
        String flag = delim+parameter+delim;
        try {
            Scanner scan = new Scanner(new File(file));
            int line = 0;
            int parseLine;
            boolean notDone = true;
            
            while(scan.hasNextLine() && notDone){
            String temp = scan.nextLine();
            line++;
                if(temp.contains(flag)){
                    return temp;
                }
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ParsingHelper.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        return "PARAMNOTFOUND";
    }

    /**
    *        parseStringsForLine parses a string list based on the parameter and returns the line 
    *   containing the parameter.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @return string parsed
    */
    public String parseStringsForLine(ArrayList<String> file,String parameter){
        String flag = delim+parameter+delim;
        String temp; 
        for(int u = 0; u < file.size();u++){        
            temp = file.get(u);
            if(temp.contains(flag)){
                return temp;
            }
        }
        return "PARAMNOTFOUND";
    }
    
    
    /**
    *        parseFileIntoStringList parses a file into a string list array.
    * 
    * @param filename - filename of file to parse
    * @return array of strings in file
    */
    public ArrayList<String> parseFileIntoStringList(String file){
        ArrayList<String> fileString = new ArrayList<String>();
        try {
            Scanner scan = new Scanner(new File(file));
            
            while(scan.hasNextLine()){
                String temp = scan.nextLine();
                fileString.add(temp);
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ParsingHelper.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        return fileString;
    }
    
    /**
    *        parse parses a file based on the parameter and returns the string of data 
    *   found and to be parsed further at input column.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return string parsed
    */
    public String parse(String input, String parameter, int locationOffset){
        String searching = delim+parameter+delim;
           
        String[] temp = input.split(searching);
        //System.out.println("Searching for: "+parameter+"    Passing: "+temp[1]);
        try{
            temp = (temp[1]).split("[ ]+");
        }catch(ArrayIndexOutOfBoundsException er){
            System.err.println("String in : "+temp[0]+"    input: "+input
                    +"   searching: "+searching);
        }
        
        ArrayList<String> arr = new ArrayList<String>();

        // clear invalid input
        for(int u = 0; u < temp.length;u++){
            // a regular number
            if(temp[u].matches("[0-9.]*") && !temp[u].equals("")){arr.add(temp[u]);}
            // a boolean
            if(temp[u].contains("true") || temp[u].contains("false")){arr.add(temp[u]);}
        }

        if(arr.size()==0){
            for(int u = 0; u < temp.length;u++){
                // a regular number
                if(!temp[u].equals("")){arr.add(temp[u]);}
            }
        }

        //System.out.println("Passing: "+temp[locationOffset]);
        return (arr.get(locationOffset));
    }
    
    /**
    *        parseIntoDblArray parses a file based on the parameter and returns the string of data 
    *   found based on parameter and passes all double data as array.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @return array parsed
    */
    public ArrayList<Double> parseIntoDblArray(String filename, String parameter){
        String searching = delim+parameter+delim;
        String input = parseFileForLine(filename,parameter);
        
        String[] temp = input.split(searching);
        //System.out.println("Searching for: "+parameter+"    Passing: "+temp[1]);
        try{temp = (temp[1]).split("[ ]+");

        }catch(ArrayIndexOutOfBoundsException er){
            System.err.println("String in : "+temp[0]+"     parameter: "+parameter
                    +"    file: "+filename);
        }
        
        ArrayList<String> arr = new ArrayList<String>();

        // clear invalid input
        for(int u = 0; u < temp.length;u++){
            // a regular number
            if(temp[u].matches("[0-9.]*") && !temp[u].equals("")){arr.add(temp[u]);}
        }

        ArrayList<Double> arrDbl = new ArrayList<Double>();
        for(int u = 0; u < arr.size();u++){
            arrDbl.add(new Double(arr.get(u)));
        }
        
        return arrDbl;
    }

        
    /**
    *        parseIntoDblArray parses a file based on the parameter and returns the string of data 
    *   found based on parameter and passes all double data as array.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @return array parsed
    */
    public ArrayList<Double> parseIntoDblArray(ArrayList<String> filename, String parameter){
        String searching = delim+parameter+delim;
        String input = parseStringsForLine(filename,parameter);
        
        String[] temp = input.split(searching);
        //System.out.println("Searching for: "+parameter+"    Passing: "+temp[1]);
        try{temp = (temp[1]).split("[ ]+");

        }catch(ArrayIndexOutOfBoundsException er){
            System.err.println("String in : "+temp[0]+"     parameter: "+parameter
                    +"    file: "+filename);
        }
        
        ArrayList<String> arr = new ArrayList<String>();

        // clear invalid input
        for(int u = 0; u < temp.length;u++){
            // a regular number
            if(temp[u].matches("[0-9.]*") && !temp[u].equals("")){arr.add(temp[u]);}
        }

        ArrayList<Double> arrDbl = new ArrayList<Double>();
        for(int u = 0; u < arr.size();u++){
            arrDbl.add(new Double(arr.get(u)));
        }
        
        return arrDbl;
    }
    
    /**
    *        parseIntoArray parses a file based on the parameter and returns the data as
    *   an array for further parsing.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @return string parsed
    */
    public ArrayList<String> parseIntoArray(String filename, String parameter){
        String searching = delim+parameter+delim;
        String input = parseFileForLine(filename,parameter);
        
        String[] temp = input.split(searching);
        if(temp.length <= 1 || temp == null){
            System.err.println("Parameter | "+searching+" |  in file : "+filename
                    +"  not parsing properly | size: "+temp.length);
        }
        //System.out.println("Searching for: "+parameter+"    Passing: "+temp[1]);
        temp = (temp[1]).split("[ ]+");

        ArrayList<String> arr = new ArrayList<String>();
        arr.add(searching);
        // clear invalid input
        for(int u = 0; u < temp.length;u++){
            Pattern p = Pattern.compile("^\\-?([0-9]{0,3}(\\,?[0-9]{3})*(\\.?[0-9]*))");   // -1,123.11
            Matcher m = p.matcher(temp[u]);
            // a regular number
            if(m.matches()){arr.add(temp[u]);}
            // a boolean
            if(temp[u].contains("true") || temp[u].contains("false")){arr.add(temp[u]);}
        }

        if(arr.size()==0){
            for(int u = 0; u < temp.length;u++){
                // a regular number
                if(!temp[u].equals("")){arr.add(temp[u]);}
            }
        }
        
        return arr;
    }

        /**
    *        parseIntoArray parses a file based on the parameter and returns the data as
    *   an array for further parsing.
    * 
    * @param filename - filename of file to parse
    * @param parameter - parameter to search for and use to parse
    * @return string parsed
    */
    public ArrayList<String> parseIntoArray(ArrayList<String> filename, String parameter){
        String searching = delim+parameter+delim;
        String input = parseStringsForLine(filename,parameter);
        
        String[] temp = input.split(searching);
        //System.out.println("Searching for: "+parameter+"    Passing: "+temp[1]);
        temp = (temp[1]).split("[ ]+");

        ArrayList<String> arr = new ArrayList<String>();
        arr.add(searching);
        // clear invalid input
        for(int u = 0; u < temp.length;u++){
            Pattern p = Pattern.compile("^\\-?([0-9]{0,3}(\\,?[0-9]{3})*(\\.?[0-9]*))");   // -1,123.11
            Matcher m = p.matcher(temp[u]);
            // a regular number
            if(m.matches()){arr.add(temp[u]);}
            // a boolean
            if(temp[u].contains("true") || temp[u].contains("false")){arr.add(temp[u]);}
        }

        if(arr.size()==0){
            for(int u = 0; u < temp.length;u++){
                // a regular number
                if(!temp[u].equals("")){arr.add(temp[u]);}
            }
        }
        
        return arr;
    }
    
    /**
    *        parseInt parses a string based on the parameter and returns integer parsed.
    * 
    * @param input - input string to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return int parsed
    */
    public int parseInt(String input, String parameter, int locationOffset){
        String param = parse(input,parameter,locationOffset);
        
        return (new Integer(param));
    }

    /**
    *        parseDouble parses a string based on the parameter and returns double parsed.
    * 
    * @param input - input string to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return double parsed
    */
    public double parseDouble(String input, String parameter, int locationOffset){
        String param = parse(input,parameter,locationOffset);
        
        return (new Double(param));
    }

    /**
    *        parseBoolean parses a string based on the parameter and returns boolean parsed.
    * 
    * @param input - input string to parse
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @return boolean parsed
    */
    public boolean parseBoolean(String input, String parameter, int locationOffset){
        String param = parse(input,parameter,locationOffset);
        if(param.contains("true")){
            return true;
        }else if(param.contains("false")){
            return false;
        }else{
            System.err.println("BOOLEAN NOT PARSED FROM: "+input);
            return false;
        }   
    }
    
    /**
    *        editFile edits a file using the parameter to parse the data and make the edit.
    * 
    * @param filename- file to edit
    * @param parameter - parameter to search for and use to parse
    * @param locationOffset - column of data to return
    * @param edit - edit to be made
    * @param maxOffset - last column to add
    */
    public void editFile(String filename, String parameter, int locationOffset, String edit, int maxOffset){
        ArrayList<String> edited = parseIntoArray(filename, parameter);
        if(edited.size() == 1){System.err.println("ERROR: "+filename+"   NOT PARSED PROPERLY FOR PARAM: "+parameter);}
        edited.set(locationOffset+1, edit);
        
        String fnameTemp = filename.substring(0, filename.length()-4)+"Temp"+".txt";
            
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fnameTemp));
            Scanner scanner = new Scanner(new File(filename));
            while(scanner.hasNextLine()) {
                String currLine = scanner.nextLine();
                if(currLine.contains(delim+parameter+delim)){
                    currLine = edited.get(0);
                    for(int u = 1;u < edited.size();u++){
                        if(u <= maxOffset){currLine =  currLine+"  "+edited.get(u);}
                    }
                }
                out.println(currLine);
            }

            out.close();

            (new File(filename)).delete();
            copyFile(fnameTemp, filename);
            (new File(fnameTemp)).delete();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    
    /**
    *      copyFile copies the original file in source file
    * 
    * @param srFile - source file
    * @param dtFile - destination file
    */    
    public void copyFile(String srFile, String dtFile){
        try{
    	  File f1 = new File(srFile);
    	  File f2 = new File(dtFile);
    	  InputStream in = new FileInputStream(f1);
    	  
    	  //For Append the file.
          //  OutputStream out = new FileOutputStream(f2,true);

    	  //For Overwrite the file.
    	  OutputStream out = new FileOutputStream(f2);

    	  byte[] buf = new byte[1024];
    	  int len;
    	  while ((len = in.read(buf)) > 0){
    	  out.write(buf, 0, len);
    	  }
    	  in.close();
    	  out.close();
    	  }
    	  catch(FileNotFoundException ex){
    	  System.out.println(ex.getMessage() + " in the specified directory.");
    	  System.exit(0);
    	  }
    	  catch(IOException e){
    	  System.out.println(e.getMessage());  
        }
    }

    /*
    *   test is the testing class.
    */
    public void test(){
        System.out.println("Result: "+parseInt(testStr,"temperature",0));
        System.out.println("Result of file: "+parseFileString(testLoc,"mcalgo", 0));
        editFile(testLoc, "field", 0, "38474.3",2);
    }
    
    public static void main(String[] args) {
        ParsingHelper parse = new ParsingHelper();
        parse.test();
    
    }
}
