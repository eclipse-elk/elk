/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.io;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.List;

/**
 * This filter can be used to filter files based on their extension. It can be for example used to get all
 * {@code .elkg} files in a directory.
 */
public class FileExtensionFilter implements FileFilter {

    /** The extension of the files that should be accepted, including the period. */
    private final List<String> extensions = new ArrayList<>();

    /**
     * Creates a new file that only accepts files with the given extension.
     * 
     * @param extension
     *            the extension to filter for, either with or without the leading period.
     */
    public FileExtensionFilter(final String extension) {
        this(new String[] {extension});
    }
    
    /**
     * Creates a new filter that only accepts files with one of the given extensions.
     * 
     * @param extensions
     *            the extension to filter for, either with or without the leading period.
     */
    public FileExtensionFilter(final String... extensions) {
        for (String ext : extensions) {
            if (ext.startsWith(".")) {
                this.extensions.add(ext);
            } else {
                this.extensions.add("." + ext);
            } 
        }
    }

    @Override
    public boolean accept(final File pathname) {
        if (!pathname.exists()) {
            return false;
        } else if (pathname.isDirectory()) {
            // Directories are accepted because the framework can look for files recursively
            return true;
        } else {
            String fileName = pathname.getName();
            return extensions.stream().anyMatch(ext -> fileName.endsWith(ext));
        }
    }

}
