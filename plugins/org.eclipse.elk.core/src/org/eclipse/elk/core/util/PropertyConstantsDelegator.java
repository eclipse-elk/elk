/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.util.HashMap;
import java.util.Map;

import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphElementAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.elk.graph.properties.Property;

/**
 * Allows to reroute access to properties through other property constants which may provide different default values.
 * This class is a workaround for the fact that algorithm-independent layout code only has {@link CoreOptions} property
 * constants to work with, along with the defaults configured therein. An algorithm can change the default value, which
 * is reflected in the algorithm's set of property constants, which (and this is the big problem) the independent code
 * doesn't know about. Instances of this class can be supplied by layout algorithms to independent code to replace its
 * {@link CoreOptions} property constants with the algorithm-specific constants.
 * 
 * <p>This class can be used in one of two ways:</p>
 * <ol>
 *   <li>Manual configuration. Call {@link #createEmpty()} to create an instance without property delegates and add
 *       delegates manually by calling {@link #addDelegate(IProperty)}.</li>
 *   <li>Automatic configuration. Call {@link #createForLayoutAlgorithmWithId(String)} to create an instance
 *       initialized with the property constants defined for the algorithm with the given ID. Requires the layout meta
 *       data service to be properly initialized.</li>
 * </ol>
 * 
 * <p><b>TODO</b> The implementation of the create methods would be much easier if layout algorithms registered their
 * property constants instead of option IDs and default values. We could then simply use the constants directly instead
 * of having to create new ones here. It would also be nice to allow a new object to be initialized from a layout
 * meta data provider, if the provider interface were extended to supply a list (or set) of supported property
 * constants. That would make the whole thing independent of the meta layout thingy and keep it working on the
 * plain Java level as well.</p>
 */
public final class PropertyConstantsDelegator {
    
    /** Set of property constants to be used for accessing their respective property. */
    private Map<IProperty<?>, IProperty<?>> propertyDelegates = new HashMap<>();

    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    /**
     * Private constructor. Use create methods instead.
     */
    private PropertyConstantsDelegator() {
        
    }
    
    /**
     * Creates and returns an empty {@link PropertyConstantsDelegator} to be configured manually.
     */
    public static PropertyConstantsDelegator createEmpty() {
        return new PropertyConstantsDelegator();
    }
    
    /**
     * Creates a new {@link PropertyConstantsDelegator} and initializes it with the layout options known to the given
     * algorithm.
     * 
     * @param algorithmData
     *            algorithm data to initialize the delegator with.
     * @return intiailzed delegator.
     */
    public static PropertyConstantsDelegator createForLayoutAlgorithmData(final LayoutAlgorithmData algorithmData) {
        PropertyConstantsDelegator delegator = new PropertyConstantsDelegator();
        
        // Delegate to property constants for all known options
        for (String optionId : algorithmData.getKnownOptionIds()) {
            delegator.addDelegate(new Property<Object>(optionId, algorithmData.getDefaultValue(optionId)));
        }
        
        return delegator;
    }
    
    /**
     * Creates a new {@link PropertyConstantsDelegator} and initializes it with the layout options known to the
     * algorithm with the given ID. If an algorithm with the ID is not registered with the
     * {@link LayoutMetaDataService}, an empty delegator is returned.
     * 
     * @param algorithmId
     *            identifier of the layout algorithm to initialize the delegator with.
     * @return intialized delegator.
     */
    public static PropertyConstantsDelegator createForLayoutAlgorithmWithId(final String algorithmId) {
        LayoutAlgorithmData algorithmData = LayoutMetaDataService.getInstance().getAlgorithmData(algorithmId);
        return algorithmData != null
                ? createForLayoutAlgorithmData(algorithmData)
                : createEmpty();
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    /**
     * Configures the given property constant to be used to access its property.
     * 
     * @param delegate
     *            the property constant to be used.
     * @return this object for method chaining.
     */
    public PropertyConstantsDelegator addDelegate(final IProperty<?> delegate) {
        propertyDelegates.put(delegate, delegate);
        return this;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters

    /**
     * Returns the delegate configured for the given property or the property itself, if there is no delegate.
     * 
     * @param property the property for which to look for a delegate.
     * @return the delegate or the property itself, if there is no delegate.
     */
    public <T> IProperty<T> getPropertyOrDelegate(final IProperty<T> property) {
        // The following cast should be safe since constants whose property IDs match should have a matching
        // type as well
        @SuppressWarnings("unchecked")
        IProperty<T> actualProperty = propertyDelegates.containsKey(property)
                ? (IProperty<T>) propertyDelegates.get(property)
                : property;
        return actualProperty;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Property Access
    
    /**
     * Calls {@link IPropertyHolder#getProperty(IProperty)} on the given property holder either using the given property
     * constant or a delegate, if one is configured.
     * 
     * @param propertyHolder
     *            the property holder to call {@link IPropertyHolder#getProperty(IProperty)} on.
     * @param property
     *            the property to access.
     * @return the property value.
     */
    public <T> T getProperty(final IPropertyHolder propertyHolder, final IProperty<T> property) {
        return propertyHolder.getProperty(getPropertyOrDelegate(property));
    }
    
    /**
     * Calls {@link IPropertyHolder#getProperty(IProperty)} on the given adapter either using the given property
     * constant or a delegate, if one is configured.
     * 
     * @param adapter
     *            the adater to call {@link IPropertyHolder#getProperty(IProperty)} on.
     * @param property
     *            the property to access.
     * @return the property value.
     */
    public <T> T getProperty(final GraphElementAdapter<?> adapter, final IProperty<T> property) {
        return adapter.getProperty(getPropertyOrDelegate(property));
    }
    
    /**
     * Calles {@link IndividualSpacings#getIndividualOrInherited(ElkNode, IProperty)} on the given node either using the
     * given property constants or a delegate, if one is configured.
     * 
     * @param node
     *            the node to retrieve the property value from.
     * @param property
     *            the property to access.
     * @return the property value.
     */
    public <T> T getIndividualOrInheritedProperty(final ElkNode node, final IProperty<T> property) {
        return IndividualSpacings.getIndividualOrInherited(node, getPropertyOrDelegate(property));
    }
    
    /**
     * Calles {@link IndividualSpacings#getIndividualOrInherited(ElkNode, IProperty)} on the given node adapter either
     * using the given property constants or a delegate, if one is configured.
     * 
     * @param nodeAdapter
     *            the node adapter to retrieve the property value from.
     * @param property
     *            the property to access.
     * @return the property value.
     */
    public <T> T getIndividualOrInheritedProperty(final NodeAdapter<?> nodeAdapter, final IProperty<T> property) {
        return IndividualSpacings.getIndividualOrInherited(nodeAdapter, getPropertyOrDelegate(property));
    }
    
}
