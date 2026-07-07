/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.io.StringWriter;
import java.util.Collections;

// elkjs-exclude-start
import org.eclipse.elk.core.util.persistence.ElkGraphResource;
// elkjs-exclude-end
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.common.util.URI;
import org.eclipse.emf.ecore.resource.Resource;
import org.eclipse.emf.ecore.resource.ResourceSet;
import org.eclipse.emf.ecore.resource.impl.ResourceSetImpl;

/**
 * A logged graph, including relevant information about the graph. The actual data type of a logged graph can be either
 * an {@link org.eclipse.elk.graph.ElkNode} (for graphs of type {@link Type#ELK}, or a string (for the other graph
 * types).
 */
public final class LoggedGraph {

    /**
     * Possible types for a logged graph.
     */
    public static enum Type {
        /** Type for elk graphs. */
        ELK(ElkNode.class, "elkg"),
        /** Type for json graphs. */
        JSON(String.class, "json"),
        /** Type for dot graphs. */
        DOT(String.class, "dot"),
        /** Type for SVG images. */
        SVG(String.class, "svg");
        
        /** Expected data type a graph of this type is expected to be delivered in. */
        private final Class<?> expectedType;
        /** Extension of files this type of graph should be stored in. */
        private final String fileExtension;
        
        Type(final Class<?> expectedType, final String fileExtension) {
            this.expectedType = expectedType;
            this.fileExtension = fileExtension;
        }
        
        /**
         * Checks if {@code o}'s type is the same or a subtype of the type we would expect.
         * 
         * @param o
         *            the object whose type compatibility to check.
         * @return {@code true} if the type is compatible.
         */
        // elkjs-exclude-start
        public boolean isTypeCompatible(final Object o) {
            return expectedType.isAssignableFrom(o.getClass());
        }
        // elkjs-exclude-end

        /**
         * Throws a {@link ClassCastException} if {@link #isTypeCompatible(Object)} returns {@code false} for the
         * object.
         * 
         * @param o
         *            the object whose type compatibility to check.
         * @throws ClassCastException if the type is unexpected.
         */
        public void checkTypeCompatibility(final Object o) {
            // elkjs-exclude-start
            if (!isTypeCompatible(o)) {
                throw new ClassCastException("Type " + o.getClass().getName() + " incompatible for " + this.name());
            }
            // elkjs-exclude-end
        }
        
        /**
         * Returns the file extension for graphs of this type.
         */
        public String getFileExtension() {
            return fileExtension;
        }
    }
    

    /** A graph. */
    private final Object graph;
    /** Type of the graph. */
    private final Type graphType;
    /** Tag for identifying the graph. */
    private final String tag;
    

    /**
     * Creates a new logged graph.
     * 
     * @throws ClassCastException
     *             if the type of {@code graph} does not conform to the type expected for {@code graphType}.
     */
    public LoggedGraph(final Object graph, final String tag, final Type graphType) {
        graphType.checkTypeCompatibility(graph);
        
        this.graph = graph;
        this.tag = tag;
        this.graphType = graphType;
    }

    /**
     * Returns the graph.
     */
    public Object getGraph() {
        return graph;
    }

    /**
     * Returns the graph type.
     */
    public Type getGraphType() {
        return graphType;
    }

    /**
     * Returns a string that will describe the graph.
     */
    public String getTag() {
        return tag;
    }
    
    /**
     * Serializes the graph to a string for it to be saved to a file.
     */
    public String serialize() {
        // Depending on our type, different things are to be done
        switch (graphType) {
        // elkjs-exclude-start
        case ELK:
            ResourceSet resourceSet = new ResourceSetImpl();
            Resource resource = resourceSet.createResource(URI.createFileURI("dummy.elkg"));
            
            // The resource should be an ElkGraphResource, but be save here...
            if (resource instanceof ElkGraphResource) {
                resource.getContents().add((ElkNode) graph);
                
                try {
                    StringWriter stringWriter = new StringWriter();
                    ((ElkGraphResource) resource).save(stringWriter, Collections.emptyMap());
                    return stringWriter.toString();
                    
                } catch (Exception e) {
                    // ignore the exception and abort the layout graph exporting
                }
            }
            
            return "Unexpected problem serializing ELK Graph.";

        // elkjs-exclude-end
        default:
            return graph.toString();
        }
    }

}
