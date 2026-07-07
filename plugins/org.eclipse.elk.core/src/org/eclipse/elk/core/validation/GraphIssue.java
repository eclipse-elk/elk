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

import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkGraphElement;

/**
 * Descriptor of an issue found while validating the structure or the configuration of a graph.
 */
public class GraphIssue {

    /** Enumeration of issue severities. */
    public static enum Severity {
        /** An error means the layout process is aborted and the error is reported to the user. */
        ERROR("Error"),
        /** A graph with warnings but no errors can still be processed. */
        WARNING("Warning");
        
        /** User-readable name of this severity. */
        private final String userString;
        
        Severity(final String userString) {
            this.userString = userString;
        }
        
        /**
         * Returns a user-friendly string for this severity.
         */
        public String getUserString() {
            return userString;
        }
    }

    /** The graph element to which the issue applies. */
    private final ElkGraphElement element;

    /** A message to be shown to users. */
    private final String message;

    /** The severity of the issue. */
    private final Severity severity;

    /**
     * Create a graph issue.
     * 
     * @param element
     *            The graph element to which the issue applies; may be {@code null}
     * @param message
     *            A message to be shown to users
     * @param severity
     *            The severity of the issue
     */
    public GraphIssue(final ElkGraphElement element, final String message, final Severity severity) {
        if (message == null || severity == null) {
            throw new NullPointerException();
        }
        this.element = element;
        this.message = message;
        this.severity = severity;
    }

    /**
     * The graph element to which the issue applies. May be {@code null}.
     */
    public ElkGraphElement getElement() {
        return element;
    }

    /**
     * A message to be shown to users.
     */
    public String getMessage() {
        return message;
    }

    /**
     * The severity of the issue.
     */
    public Severity getSeverity() {
        return severity;
    }
    
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof GraphIssue) {
            GraphIssue other = (GraphIssue) obj;
            return this.element == other.element
                    && this.message.equals(other.message)
                    && this.severity == other.severity;
        }
        return false;
    }
    
    @Override
    public int hashCode() {
        return element.hashCode() ^ message.hashCode() ^ severity.hashCode();
    }
    
    @Override
    public String toString() {
        StringBuilder result = new StringBuilder();
        result.append(severity).append(": ").append(message).append(" (");
        ElkUtil.printElementPath(element, result);
        result.append(")");
        return result.toString();
    }

}
