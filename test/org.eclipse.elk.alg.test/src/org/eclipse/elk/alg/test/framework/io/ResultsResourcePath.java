/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
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

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.google.common.base.Strings;

/**
 * The location of a results file relative to the path specified by the {@code RESULTS_PATH} system property.
 */
public class ResultsResourcePath extends AbstractPropertyDependentResourcePath {
    
    /** The system property or environment variable that contains our base path. */
    public static final String PATH_PROPERTY = "RESULTS_PATH";

    /**
     * Creates a new instance that points to the given resource.
     * 
     * @param relativePath
     *            path to the resulting file, relative to the {@code RESULTS_PATH} system property.
     */
    public ResultsResourcePath(final String relativePath) {
        if (Strings.isNullOrEmpty(relativePath)) {
            throw new IllegalArgumentException("The file path cannot be empty.");
        }
        
        initialize(basePathForProperty(PATH_PROPERTY), relativePath);
    }
    
    /**
     * Ensures that the results directory actually exists.
     * 
     * @throws IOException if the path couldn't be created.
     */
    public static void ensureResultsPathExists() throws IOException {
        String path = basePathForProperty(PATH_PROPERTY);
        Files.createDirectories(Paths.get(path));
    }

}
