/**
 * Adapated from examples provided by NCBI
 *
 * // original @author  Jim Mork
 * @author Alice Gutteridge
 * // original @version 1.0, September 18, 2006
**/

import gov.nih.nlm.nls.skr.*;
import java.io.*;

public class MetaMapCaller {

	public static void main(String[] args) throws InterruptedException {
		// username must be args[1], password must be args[2]
		String username = args[1];
		String password = args[2];
		GenericObject myGenericObj = new GenericObject(username, password);

		// filename path must be args[0], valid email address must be args[3]
		myGenericObj.setFileField("UpLoad_File", args[0]);
		myGenericObj.setField("Email_Address", args[3]);
		myGenericObj.setField("Batch_Command", "MTI -opt1L_DCMS -E");
		myGenericObj.setField("BatchNotes", "SKR Web API test");
		myGenericObj.setField("SilentEmail", true);

		try {
			String results = myGenericObj.handleSubmission();
			System.out.println(results);
			System.out.println("END");
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
		}
	}
}