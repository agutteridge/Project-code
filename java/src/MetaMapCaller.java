/**
 * Adapated from examples provided by NCBI
 * Example program for submitting a new Generic Batch with Validation job
 * request to the Scheduler to run. You will be prompted for your username and
 * password and if they are alright, the job is submitted to the Scheduler and
 * the results are returned in the String "results" below.
 *
 * This example shows how to setup a basic Generic Batch with Validation job
 * with a small file (sample.txt) with ASCII MEDLINE formatted citations as
 * input data. You must set the Email_Address variable and use the UpLoad_File
 * to specify the data to be processed.  This example also shows the user
 * setting the SilentEmail option which tells the Scheduler to NOT send email
 * upon completing the job.
 *
 * This example is set to run the MTI (Medical Text Indexer) program using
 * the -opt1L_DCMS and -E options. You can also setup any environment variables
 * that will be needed by the program by setting the Batch_Env field.
 * The "-E" option is required for all of the various SKR tools (MetaMap,
 * SemRep, and MTI), so please make sure to add the option to your command!
 * 
 * // original @author  Jim Mork
 * @author Alice Gutteridge
 * // original @version 1.0, September 18, 2006
**/

import gov.nih.nlm.nls.skr.*;
import java.io.*;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Scanner;
import org.json.simple.JSONValue;

public class MetaMapCaller {

	public static void main(String[] args) throws InterruptedException {
		GenericObject myGenericObj = new GenericObject();

		// filename must be args[0], valid email address must be args[1]
		myGenericObj.setFileField("UpLoad_File", args[0]);
		myGenericObj.setField("Email_Address", args[1]);
		myGenericObj.setField("Batch_Command", "MTI -opt1L_DCMS -E");
		myGenericObj.setField("BatchNotes", "SKR Web API test");
		myGenericObj.setField("SilentEmail", true);

		try {
			String results = myGenericObj.handleSubmission();

			String[] resultArray = results.split("\n");
			LinkedList<String> list = new LinkedList<String>(Arrays.asList(resultArray));
			StringWriter out = new StringWriter();
			// using JSON.simple library 
			JSONValue.writeJSONString(list, out);
			String jsonText = out.toString();
			System.out.print(jsonText);

		} catch (RuntimeException ex) {
			System.err.println("");
			System.err.print("An ERROR has occurred while processing your");
			System.err.println(" request, please review any");
			System.err.print("lines beginning with \"Error:\" above and the");
			System.err.println(" trace below for indications of");
			System.err.println("what may have gone wrong.");
			System.err.println("");
			System.err.println("Trace:");
			ex.printStackTrace();
		} catch (IOException io) {
			System.err.println("An exception was thrown whilst writing JSON.");
			System.err.println("The input is incompatible with the JSONValue.writeJSONString operation.");
		}
	}
}