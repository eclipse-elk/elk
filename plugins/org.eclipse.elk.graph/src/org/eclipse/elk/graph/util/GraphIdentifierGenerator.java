/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.graph.util;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;

import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.emf.ecore.EObject;

import com.google.common.base.Strings;
import com.google.common.collect.Iterators;
import com.google.common.collect.Sets;

/**
 * Generates identifiers for graph elements where missing. Inside ELK, this class is mainly used to generate
 * proper identifiers before a graph is serialized using the textual graph languages, which uses identifiers for
 * elements referencing each other.
 */
public final class GraphIdentifierGenerator {
    
    private ElkNode graph;
    private boolean validate = false;
    private boolean generate = false;
    private boolean unique = false;
    
    private Random random = null;
    
    /**
     * @param graph the graph for which to generate identifiers. 
     * @return a new {@link GraphIdentifierGenerator} instance that can be configured and finally executed. 
     * 
     * @see #assertExists()
     * @see #assertValid()
     * @see #assertUnique()
     * @see #execute()
     */
    public static GraphIdentifierGenerator forGraph(final ElkNode graph) {
        return new GraphIdentifierGenerator(graph);
    }
    
    /**
     * Ensures that existing identifiers are valid instances of identifiers as defined by the ELKT format.
     */
    public GraphIdentifierGenerator assertValid() {
        validate = true;
        return this;
    }
    
    /**
     * Generates identifiers for the graph's elements, where missing.
     */
    public GraphIdentifierGenerator assertExists() {
        generate = true;
        return this;
    }
    
    /**
     * Makes sure existing identifiers are unique by appending a "_" to non-unique ones. 
     */
    public GraphIdentifierGenerator assertUnique() {
        unique = true;
        return this;
    }
    
    /**
     * Execute this generator. 
     */
    public void execute() {
        if (validate) {
            validateIdentifiers(graph);
        }
        
        if (generate) {
            generateIdentifiers(graph);
        }
        
        if (unique) {
            assertAllIdsUnique(graph);
        }
    }
    
    /**
     * Enumeration of possible graph elements that can receive identifiers.
     */
    private static enum ElementType {
        NODE("N"),
        PORT("P"),
        EDGE("E"),
        EDGE_SECTION("ES"),
        LABEL("L");
        
        private String elementPrefix;
        
        ElementType(final String prefix) {
            elementPrefix = prefix;
        }
    }
    
    /** Array of highest identifiers generated so far for each kind of element. */
    private int[] currentIDs = new int[ElementType.values().length];
    /** Set of existing identifiers in the graph for collision detection. */
    private Set<String> existingIdentifiers = new HashSet<>();
    
    
    /**
     * Constructor is only called from the inside.
     */
    private GraphIdentifierGenerator(final ElkNode graph) {
        this.graph = graph;
        Iterator<EObject> iterator = graph.eAllContents();
        while (iterator.hasNext()) {
            EObject currentEObject = iterator.next();
            
            if (currentEObject instanceof ElkGraphElement) {
                ElkGraphElement element = (ElkGraphElement) currentEObject;
                if (element.getIdentifier() != null && !element.getIdentifier().trim().isEmpty()) {
                    existingIdentifiers.add(element.getIdentifier());
                }
            } else if (currentEObject instanceof ElkEdgeSection) {
                ElkEdgeSection section = (ElkEdgeSection) currentEObject;
                if (section.getIdentifier() != null && !section.getIdentifier().trim().isEmpty()) {
                    existingIdentifiers.add(section.getIdentifier());
                }
            }
        }
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Identifier Validation
    
    /**
     * Ensures valid identifiers in the object tree rooted at the given element.
     */
    private void validateIdentifiers(final EObject element) {
        new ElkGraphSwitch<Boolean>() {
            
            @Override
            public Boolean caseElkNode(final ElkNode node) {
                validateIdentifier(node);
                
                node.getLabels().stream().forEach(l -> validateIdentifiers(l));
                node.getPorts().stream().forEach(p -> validateIdentifiers(p));
                node.getContainedEdges().stream().forEach(e -> validateIdentifiers(e));
                node.getChildren().stream().forEach(c -> validateIdentifiers(c));
                return true;
            }
            
            @Override
            public Boolean caseElkPort(final ElkPort port) {
                validateIdentifier(port);
                
                port.getLabels().stream().forEach(l -> validateIdentifiers(l));
                return true;
            }
            
            @Override
            public Boolean caseElkLabel(final ElkLabel label) {
                validateIdentifier(label);
                
                label.getLabels().stream().forEach(l -> validateIdentifiers(l));
                return true;
            }
            
            @Override
            public Boolean caseElkEdge(final ElkEdge edge) {
                validateIdentifier(edge);
                
                edge.getLabels().stream().forEach(l -> validateIdentifiers(l));
                edge.getSections().stream().forEach(s -> validateIdentifiers(s));
                return true;
            }
            
            @Override
            public Boolean caseElkEdgeSection(final ElkEdgeSection section) {
                validateIdentifier(section);
                return true;
            }
            
        }.doSwitch(element);
    }
    
    private void validateIdentifier(final ElkGraphElement element) {
        String validIdentifier = validateIdentifier(element.getIdentifier());
        if (validIdentifier != null) {
            element.setIdentifier(validIdentifier);
        }
    }
    
    private void validateIdentifier(final ElkEdgeSection section) {
        String validIdentifier = validateIdentifier(section.getIdentifier());
        if (validIdentifier != null) {
            section.setIdentifier(validIdentifier);
        }
    }
    
    /**
     * Returns a validated version of the given identifier or {@code null} if the identifier was already valid (or
     * empty or {@code null} itself).
     */
    private String validateIdentifier(final String identifier) {
        if (Strings.isNullOrEmpty(identifier)) {
            return null;
        }
        
        boolean valid = true;
        char[] chars = identifier.toCharArray();
        
        for (int i = 0; i < chars.length; i++) {
            // Rule to validate: ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( 'a' .. 'z' | 'A' .. 'Z' | '_' | '0' .. '9' )*
            if (!(chars[i] >= 'A' && chars[i] <= 'Z')
                    && !(chars[i] >= 'a' && chars[i] <= 'z')
                    && !(chars[i] == '_')
                    && !(i > 0 && chars[i] >= '0' && chars[i] <= '9')) {
                
                // We have an invalid character, replace with _
                chars[i] = '_';
                valid = false;
            }
        }
        
        return valid ? null : new String(chars);
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Identifier Generation
   
    /**
     * Recursively generates identifiers for the given element and its child elements.
     * 
     * @param element the element to generate an identifier for.
     */
    private GraphIdentifierGenerator generateIdentifiers(final EObject element) {
        new ElkGraphSwitch<Boolean>() {
            
            @Override
            public Boolean caseElkNode(final ElkNode node) {
                if (node.getParent() == null) {
                    if (node.getIdentifier() == null || node.getIdentifier().trim().isEmpty()) {
                        node.setIdentifier("G1");
                    }
                } else {
                    setIdentifierIfMissing(node, ElementType.NODE);
                }
                
                node.getLabels().stream().forEach(l -> generateIdentifiers(l));
                node.getPorts().stream().forEach(p -> generateIdentifiers(p));
                node.getContainedEdges().stream().forEach(e -> generateIdentifiers(e));
                node.getChildren().stream().forEach(c -> generateIdentifiers(c));
                return true;
            }
            
            @Override
            public Boolean caseElkPort(final ElkPort port) {
                setIdentifierIfMissing(port, ElementType.PORT);
                
                port.getLabels().stream().forEach(l -> generateIdentifiers(l));
                return true;
            }
            
            @Override
            public Boolean caseElkLabel(final ElkLabel label) {
                setIdentifierIfMissing(label, ElementType.LABEL);
                
                label.getLabels().stream().forEach(l -> generateIdentifiers(l));
                return true;
            }
            
            @Override
            public Boolean caseElkEdge(final ElkEdge edge) {
                setIdentifierIfMissing(edge, ElementType.EDGE);
                
                edge.getLabels().stream().forEach(l -> generateIdentifiers(l));
                edge.getSections().stream().forEach(s -> generateIdentifiers(s));
                return true;
            }
            
            @Override
            public Boolean caseElkEdgeSection(final ElkEdgeSection section) {
                setIdentifierIfMissing(section);
                return true;
            }
            
        }.doSwitch(element);
        
        return this;
    }
    
    /**
     * Generates and sets a new identifier for the given element of the given type.
     */
    private void setIdentifierIfMissing(final ElkGraphElement element, final ElementType elementType) {
        if (element.getIdentifier() == null || element.getIdentifier().trim().isEmpty()) {
            String id = nextIdentifier(elementType);
            element.setIdentifier(id);
        }
    }
    
    /**
     * Generates and sets a new identifier for the given edge section.
     */
    private void setIdentifierIfMissing(final ElkEdgeSection section) {
        if (section.getIdentifier() == null || section.getIdentifier().trim().isEmpty()) {
            String id = nextIdentifier(ElementType.EDGE_SECTION);
            section.setIdentifier(id);
        }
    }
    
    /**
     * Returns the next unused identifier for an element of the given type.
     */
    private String nextIdentifier(final ElementType elementType) {
        String identifier;
        do {
            identifier = elementType.elementPrefix + ++currentIDs[elementType.ordinal()];
        } while (existingIdentifiers.contains(identifier));
        
        return identifier;
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Uniqueness
    
    private GraphIdentifierGenerator assertAllIdsUnique(final EObject element) {
        Set<String> knownIds = Sets.newHashSet();
        Iterator<ElkGraphElement> elementIt = Iterators.filter(element.eAllContents(), ElkGraphElement.class);
        while (elementIt.hasNext()) {
            ElkGraphElement e = elementIt.next();
            while (knownIds.contains(e.getIdentifier())) {
                String newId = e.getIdentifier() + "_g" + fourDigitPaddedRandomNumber();
                e.setIdentifier(newId);
            }
            knownIds.add(e.getIdentifier());
        }
        
        return this;
    }
    
    private String fourDigitPaddedRandomNumber() {
        if (random == null) {
            random = new Random();
        }
        // SUPPRESS CHECKSTYLE NEXT 2 MagicNumber
        int rand = random.nextInt(10000);
        return padZeroes(Integer.toString(rand), 4);
    }
    
    private String padZeroes(final String s, final int length) {
        if (s.length() == length) {
            return s;
        }
        StringBuffer r = new StringBuffer(length);
        for (int i = 0; i < length - s.length(); ++i) {
            r.append("0");
        }
        r.append(s);
        return r.toString();
    }
    
}
