/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
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
package org.eclipse.elk.core.service;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Set;

import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutAlgorithmResolver;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.properties.GraphFeature;
import org.eclipse.elk.graph.properties.IProperty;

import com.google.common.collect.Iterators;
import com.google.inject.Inject;
import com.google.inject.Singleton;

/**
 * Manages layout configurations derived from implementations of {@link ILayoutConfigurationStore}.
 * 
 * <p>Subclasses of this class can be bound in an {@link ILayoutSetup} injector for customization.
 * Note that this class is marked as {@link Singleton}, which means that exactly one instance is
 * created for each injector, i.e. for each registered {@link ILayoutSetup}.</p>
 */
@Singleton
public class LayoutConfigurationManager {
    
    /**
     * Optional provider for layout configuration stores.
     */
    @Inject(optional = true)
    private ILayoutConfigurationStore.Provider configProvider;
    
    /** The graph layout engine defines the default layout algorithm to be used. */
    @Inject
    private LayoutAlgorithmResolver layoutAlgorithmResolver;
    
    /**
     * Return a layout algorithm data instance from the given configuration store. If no specific
     * algorithm is configured, a default algorithm is selected.
     * 
     * @return a layout algorithm, or {@code null} if none was found for the given configuration
     */
    public LayoutAlgorithmData getAlgorithm(final ILayoutConfigurationStore config) {
        return LayoutMetaDataService.getInstance().getAlgorithmDataBySuffixOrDefault(
                (String) getRawOptionValue(CoreOptions.ALGORITHM, config),
                layoutAlgorithmResolver.getDefaultLayoutAlgorithmID());
    }
    
    /**
     * Determine a list of options supported by the layout algorithms from the given configuration store.
     * 
     * @param config a layout configuration store associated to a graph element
     * @return a list of supported options
     */
    public Set<LayoutOptionData> getSupportedOptions(final ILayoutConfigurationStore config) {
        LayoutMetaDataService layoutDataService = LayoutMetaDataService.getInstance();
        
        Set<LayoutOptionData> optionData = new LinkedHashSet<LayoutOptionData>();
        Set<LayoutOptionData.Target> optionTargets = config.getOptionTargets();
        if (optionTargets.contains(LayoutOptionData.Target.PARENTS)) {
            LayoutAlgorithmData algoData = getAlgorithm(config);
            if (algoData != null) {
                optionData.addAll(layoutDataService.getOptionData(algoData, LayoutOptionData.Target.PARENTS));
            }
        }
        ILayoutConfigurationStore parentConfig = config.getParent();
        if (parentConfig != null) {
            LayoutAlgorithmData algoData = getAlgorithm(parentConfig);
            if (algoData != null) {
                for (LayoutOptionData.Target target : optionTargets) {
                    if (target != LayoutOptionData.Target.PARENTS) {
                        optionData.addAll(layoutDataService.getOptionData(algoData, target));
                    }
                }
            }
        }
        return optionData;
    }
    
    /**
     * Retrieve the currently valid value for the given layout option. If no value is found in the
     * given configuration store, a default value is taken.
     * 
     * @param optionData a layout option descriptor
     * @param config a layout configuration store from which to fetch values
     */
    public Object getOptionValue(final LayoutOptionData optionData, final ILayoutConfigurationStore config) {
        Object result = getRawOptionValue(optionData, config);
        if (optionData.equals(CoreOptions.ALGORITHM)) {
            LayoutAlgorithmData algoData = LayoutMetaDataService.getInstance().getAlgorithmDataBySuffixOrDefault(
                    (String) result, layoutAlgorithmResolver.getDefaultLayoutAlgorithmID());
            if (algoData != null) {
                return algoData.getId();
            }
        } else if (result != null) {
            return result;
        }
        
        if (optionData.getTargets().contains(LayoutOptionData.Target.PARENTS)) {
            // Check default value of the content layout algorithm
            LayoutAlgorithmData algoData = getAlgorithm(config);
            if (algoData != null) {
                result = algoData.getDefaultValue(optionData);
                if (result != null) {
                    return result;
                }
            }
        } else {
            // Check default value of the container layout algorithm
            ILayoutConfigurationStore parentConfig = config.getParent();
            if (parentConfig != null) {
                LayoutAlgorithmData algoData = getAlgorithm(parentConfig);
                if (algoData != null) {
                    result = algoData.getDefaultValue(optionData);
                    if (result != null) {
                        return result;
                    }
                }
            }
        }
        
        // Fall back to default value of the option itself
        result = optionData.getDefault();
        if (result != null) {
            return result;
        }
        
        // Fall back to default-default value
        return optionData.getDefaultDefault();
    }
    
    /**
     * Retrieve the currently valid value for the given layout option. Default values are not considered.
     * 
     * @param option a layout option
     * @param config a layout configuration store from which to fetch values
     */
    protected Object getRawOptionValue(final IProperty<?> option, final ILayoutConfigurationStore config) {
        Object result = config.getOptionValue(option.getId());
        
        if (result instanceof String && option instanceof LayoutOptionData) {
            return ((LayoutOptionData) option).parseValue((String) result);
        } else {
            return result;
        }
    }
    
    /**
     * Determine the parent configuration store that is responsible for the given child configuration.
     * The highest parent for which full-hierarchy layout is enabled is responsible, if it exists,
     * and otherwise it's the direct parent (which may be {@code null}).
     */
    protected ILayoutConfigurationStore getResponsibleParent(final ILayoutConfigurationStore config) {
        ILayoutConfigurationStore result = config.getParent();
        ILayoutConfigurationStore current = result;
        ILayoutConfigurationStore parent;
        while (current != null) {
            parent = current.getParent();
            if (parent != null && isFullHierarchyLayout(parent)) {
                result = parent;
            }
            current = parent;
        }
        return result;
    }
    
    /**
     * Determine whether the given configuration store is set up to do a full-hierarchy layout,
     * i.e. to process all contained hierarchy levels with a single layout algorithm invocation.
     */
    protected boolean isFullHierarchyLayout(final ILayoutConfigurationStore config) {
        Object layoutHierarchyVal = config.getOptionValue(CoreOptions.HIERARCHY_HANDLING.getId());
        HierarchyHandling layoutHierarchy = layoutHierarchyVal instanceof HierarchyHandling
                ? (HierarchyHandling) layoutHierarchyVal
                : HierarchyHandling.valueOf((String) layoutHierarchyVal);
                
        if (layoutHierarchy == HierarchyHandling.INCLUDE_CHILDREN) {
            LayoutAlgorithmData algoData = getAlgorithm(config);
            return algoData != null
                    && (algoData.supportsFeature(GraphFeature.COMPOUND)
                    || algoData.supportsFeature(GraphFeature.CLUSTERS));
        }
        return false;
    }
    
    /**
     * Return the root configuration store for the given one by going up the parent hierarchy.
     */
    protected ILayoutConfigurationStore getRoot(final ILayoutConfigurationStore config) {
        ILayoutConfigurationStore current = config;
        ILayoutConfigurationStore parent;
        do {
            parent = current.getParent();
            if (parent == null) {
                return current;
            }
            current = parent;
        } while (current != null);
        return null;
    }
    
    /**
     * Clear all non-recursive layout option values for the given configuration store.
     */
    public void clearOptionValues(final ILayoutConfigurationStore config) {
        for (String optionId : config.getAffectedOptions()) {
            config.setOptionValue(optionId, null);
        }
    }
    
    /**
     * Create a layout configurator and initialize it with the option values from the given layout
     * configuration store provider.
     */
    public <T> LayoutConfigurator createConfigurator(final LayoutMapping layoutMapping) {
        LayoutConfigurator result = new LayoutConfigurator();
        if (configProvider != null) {
            configureElement(layoutMapping.getLayoutGraph(), layoutMapping, result);
            Iterator<ElkGraphElement> allElements = Iterators.filter(layoutMapping.getLayoutGraph().eAllContents(),
                    ElkGraphElement.class);
            while (allElements.hasNext()) {
                ElkGraphElement element = allElements.next();
                configureElement(element, layoutMapping, result);
            }
        }
        return result;
    }
    
    /**
     * Transfer option values from a configuration store of a graph element to a configurator. The
     * configuration store is obtained for the diagram part that corresponds to the given graph element
     * in the layout mapping. 
     */
    protected void configureElement(final ElkGraphElement element, final LayoutMapping layoutMapping,
            final LayoutConfigurator configurator) {
        Object diagramPart = layoutMapping.getGraphMap().get(element);
        ILayoutConfigurationStore configurationStore = configProvider.get(layoutMapping.getWorkbenchPart(),
                diagramPart);
        if (configurationStore != null) {
            configureElement(element, configurationStore, configurator);
        }
    }
    
    /**
     * Transfer option values from the given configuration store to a configurator using the given graph
     * element as key.
     */
    protected void configureElement(final ElkGraphElement element, final ILayoutConfigurationStore configStore,
            final LayoutConfigurator configurator) {
        LayoutMetaDataService layoutDataService = LayoutMetaDataService.getInstance();
        for (String optionId : configStore.getAffectedOptions()) {
            Object value = configStore.getOptionValue(optionId);
            LayoutOptionData optionData = layoutDataService.getOptionData(optionId);
            if (optionData != null && value != null) {
                if (value instanceof String) {
                    value = optionData.parseValue((String) value);
                }
                configurator.configure(element).setProperty(optionData, value);
            }
        }
    }

}
