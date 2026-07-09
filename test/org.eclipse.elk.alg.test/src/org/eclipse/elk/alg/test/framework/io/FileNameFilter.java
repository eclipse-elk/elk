/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.io;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * This filter can be used to filter files based on their file name, matched against regular expressions.
 */
public class FileNameFilter implements FileFilter {

    /** Patterns to match the file names against. */
    private final List<Pattern> patterns = new ArrayList<>();

    /**
     * Creates a new file that only accepts files that match the given expression.
     * 
     * @param expression
     *            the expression to filter for.
     */
    public FileNameFilter(final String expression) {
        this(new String[] {expression});
    }
    
    /**
     * Creates a new filter that only accepts files that match at least one of the given expressions.
     * 
     * @param expressions
     *            the extension to filter for, either with or without the leading period.
     */
    public FileNameFilter(final String... expressions) {
        for (String exp : expressions) {
             patterns.add(Pattern.compile(exp));
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
            return patterns.stream().anyMatch(pattern -> pattern.matcher(fileName).matches());
        }
    }

}
