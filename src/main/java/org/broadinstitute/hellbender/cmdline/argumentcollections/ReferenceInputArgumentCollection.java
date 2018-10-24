package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.Serializable;
import java.nio.file.Path;

/**
 * An abstract ArgumentCollection for specifying a reference sequence file
 */
public abstract class ReferenceInputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Get the name of the reference file specified at the command line.
     */
    // TODO: the name (and return type) of this should change, since it used to return a String suitable
    // for use with a File constructor; it currently returns a URI string that contains a scheme, which
    // can be used to create a Path, but not a java.nio.File constructor; and it should probably just return
    // GATKInputPathSpecifier
    public abstract String getReferenceFileName();

    /**
     * Get the Path to the reference, may be null
     */
    // TODO: Not sure we need to keep this; if getReferenceFileName returned a GATKInputPathSpecifier
    // then callers can just call toPath on that just-in-time
    public Path getReferencePath() {
        return getReferenceFileName() != null ? IOUtils.getPath(getReferenceFileName()) : null;
    }
}
