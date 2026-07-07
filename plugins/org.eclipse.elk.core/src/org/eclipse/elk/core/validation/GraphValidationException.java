/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.validation;

import java.util.Collection;

/**
 * Thrown when a graph validator found at least one error.
 */
public class GraphValidationException extends RuntimeException {

    /** the serial version UID. */
    private static final long serialVersionUID = 669762537737088914L;
    
    /** Issues attached to the exception. */
    private Collection<GraphIssue> issues;
    
    /**
     * Create a graph validation exception with no parameters.
     */
    public GraphValidationException() {
        super();
    }
    
    /**
     * Create a graph validation exception with a message.
     */
    public GraphValidationException(final String message) {
        super(message);
    }
    
    /**
     * Create a graph validation exception with a message and a collection of graph issues.
     */
    public GraphValidationException(final String message, final Collection<GraphIssue> issues) {
        super(message);
        this.issues = issues;
    }
    
    /**
     * Issues attached to the exception.
     */
    public Collection<GraphIssue> getIssues() {
        return issues;
    }

}
