/*******************************************************************************
 * Copyright (c) 2016, 2020 Kiel University and others.
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
package org.eclipse.elk.core.util;

import java.util.stream.Collectors;

import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

/**
 * A property holder which can be used to define individual spacing overrides. The overrides are applied to the
 * element this object is set on, which is done through the {@link CoreOptions#SPACING_INDIVIDUAL} property.
 */
public class IndividualSpacings extends MapPropertyHolder implements IDataObject {

    /** Serialization identifier. */
    private static final long serialVersionUID = 737614242607924309L;
    
    /**
     * Constructs an empty spacings container.
     */
    public IndividualSpacings() { }
    
    /**
     * Constructs a new spacings container and copies all properties of {@code other} to the new container.
     */
    public IndividualSpacings(final IndividualSpacings other) {
        getAllProperties().putAll(other.getAllProperties());
    }
    
    /**
     * Returns the value of the given property as it applies to the given node. First checks whether an individual
     * override is set on the node that has the given property configured. If so, the configured value is returned.
     * Otherwise, the node's parent node, if any, is queried.
     * 
     * @param node
     *            the node whose property value to return.
     * @param property
     *            the property.
     * @return the individual override for the property or the default value inherited by the parent node.
     */
    public static <T> T getIndividualOrInherited(final ElkNode node, final IProperty<T> property) {
        T result = null;
        
        if (node.hasProperty(CoreOptions.SPACING_INDIVIDUAL)) {
            IPropertyHolder individualSpacings = node.getProperty(CoreOptions.SPACING_INDIVIDUAL);
            if (individualSpacings.hasProperty(property)) {
                result = individualSpacings.getProperty(property);
            }
        }
        
        // use the common value
        if (result == null && node.getParent() != null) {
            result = node.getParent().getProperty(property);
        }
        
        return result;
    }
    
    /**
     * Returns the value of the given property as it applies to the given node. First checks whether an individual
     * override is set on the node that has the given property configured. If so, the configured value is returned.
     * Otherwise, the graph the node is part of, if any, is queried.
     * 
     * @param node
     *            the node whose property value to return.
     * @param property
     *            the property.
     * @return the individual override for the property or the default value inherited by the parent node.
     */
    public static <T> T getIndividualOrInherited(final NodeAdapter<?> node, final IProperty<T> property) {
        T result = null;
        
        if (node.hasProperty(CoreOptions.SPACING_INDIVIDUAL)) {
            IPropertyHolder individualSpacings = node.getProperty(CoreOptions.SPACING_INDIVIDUAL);
            if (individualSpacings.hasProperty(property)) {
                result = individualSpacings.getProperty(property);
            }
        }
        
        // use the common value
        if (result == null && node.getGraph() != null) {
            result = node.getGraph().getProperty(property);
        }
        
        // if the result is still null, we need the property's default value
        if (result == null) {
            result = property.getDefault();
        }
        
        return result;
    }
    
    /** A (hopefully) unique separator that allows single occurrences of commas, colons, and semi-colons in-between. */
    private static final String SERIALIZED_OPTION_SEPARATOR = ";,;";
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        final String serialized = getAllProperties().entrySet().stream()
                .map(entry -> entry.getKey().getId() + ":" + entry.getValue().toString())
                .collect(Collectors.joining(SERIALIZED_OPTION_SEPARATOR));
        return serialized;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void parse(final String string) {
        if (string == null || string.isEmpty()) {
            // Nothing to do.
            return;
        }
        try {
            String[] options = string.split(SERIALIZED_OPTION_SEPARATOR);
            for (String optionString : options) {
                String[] option = optionString.split("\\:");
                LayoutOptionData optionData = LayoutMetaDataService.getInstance().getOptionDataBySuffix(option[0]);
                if (optionData == null) {
                    throw new IllegalArgumentException("Invalid option id: " + option[0]);
                }
                Object value = optionData.parseValue(option[1]);
                if (value == null) {
                    throw new IllegalArgumentException("Invalid option value: " + option[1]);
                }
                setProperty(optionData, value);
            }
        } catch (Exception e) {
            throw new IllegalArgumentException(
                    "The given string does not match the expected format for individual spacings.", e);
        }
    }

}
