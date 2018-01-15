package PfTools;
public class Profile {
    static {
        try {
            //String OS = System.getProperty("os.name").toLowerCase();
            System.loadLibrary("prf");
        } catch (Exception e) {
            System.exit(-1);
        }
    }

    protected long profilePointer;
    protected long workPointer;

    public Profile(String profile) {
			profilePointer = 0;
			workPointer = 0;
			if (!loadProfile(profile)) {
				throw new RuntimeException("Error loading profile " + profile);
			}
    }
        
    @Override
		public void finalize() {
			freeProfile();
		}

		native public boolean createJSON(String FileName);
    native private boolean loadProfile(String profileFile);
		native private void freeProfile();
}
