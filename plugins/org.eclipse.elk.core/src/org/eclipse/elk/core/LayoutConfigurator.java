/*******************************************************************************
 * Copyright (c) 2015, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core;

import java.util.List;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.util.IGraphElementVisitor;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.elk.graph.properties.MapPropertyHolder;
import org.eclipse.elk.graph.properties.Property;
import org.eclipse.elk.graph.util.ElkReflect;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * A layout configurator is a graph element visitor that applies layout option values. It can be
 * used to modify the layout configuration of a graph after it has been created, e.g. in order to
 * apply multiple layouts with different configurations. Create an instance and then use one of the
 * {@code configure(..)} methods to obtain a property holder that can be filled with values for
 * layout options.
 */
public class LayoutConfigurator implements IGraphElementVisitor {
    
    /**
     * Property attached to the shape layout of the top-level {@link ElkNode} holding a reference
     * to an additional layout configurator that shall be applied for the remaining iterations.
     * This is applicable only when multiple layout iterations are performed, e.g. for applying
     * different layout algorithms one after another. If only a single iteration is performed or
     * the current iteration is the last one, the value of this property is ignored. Otherwise
     * the referenced configurator will be added to all iterations that follow the current one.
     */
    public static final IProperty<LayoutConfigurator> ADD_LAYOUT_CONFIG =
            new Property<LayoutConfigurator>("org.eclipse.elk.addLayoutConfig");
    
    private final Map<ElkGraphElement, MapPropertyHolder> elementOptionMap = Maps.newHashMap();
    private final Map<Class<? extends ElkGraphElement>, MapPropertyHolder> classOptionMap = Maps.newHashMap();
    private boolean clearLayout = false;
    private final List<IOptionFilter> optionFilters = Lists.newArrayList();
    
    /**
     * Functional interface that allows to specify whether a certain property should be set for a certain graph element.
     */
    @FunctionalInterface
    public interface IOptionFilter {
        /**
         * @return {@code true} if and only if the passed graph element may be configured with the passed property.
         */
        boolean accept(ElkGraphElement e, IProperty<?> p);
    }
    
    /**
     * Functional interface that allows to specify whether a certain property should be set for a certain property
     * holder.
     */
    @FunctionalInterface
    public interface IPropertyHolderOptionFilter {
        /**
         * @return {@code true} if and only if the passed property holder may be configured with the passed property.
         */
        boolean accept(IPropertyHolder holder, IProperty<?> property);
    }

    /**
     * Generic filter that prevents the {@link LayoutConfigurator} from overwriting layout options that are 
     * already set for a graph element.
     */
    public static final IPropertyHolderOptionFilter NO_OVERWRITE_HOLDER = (e, p) -> !e.hasProperty(p);
    
    /**
     * Generic filter that prevents the {@link LayoutConfigurator} from overwriting layout options that are 
     * already set for a graph element.
     */
    public static final IOptionFilter NO_OVERWRITE = (e, p) -> !e.hasProperty(p);
    
    /**
     * A filter for that checks for each option whether its configured targets match the input element. Uses the
     * {@link LayoutOptionData} obtained from the {@link LayoutMetaDataService} to do so.
     */
    public static final IOptionFilter OPTION_TARGET_FILTER = (e, property) -> {
        LayoutOptionData optionData = LayoutMetaDataService.getInstance().getOptionData(property.getId());
        if (optionData != null) {
            Set<LayoutOptionData.Target> targets = optionData.getTargets();
            if (e instanceof ElkNode) {
                if (!((ElkNode) e).isHierarchical()) {
                    return targets.contains(LayoutOptionData.Target.NODES);
                } else {
                    return targets.contains(LayoutOptionData.Target.NODES)
                            || targets.contains(LayoutOptionData.Target.PARENTS);
                }
            } else if (e instanceof ElkEdge) {
                return targets.contains(LayoutOptionData.Target.EDGES);
            } else if (e instanceof ElkPort) {
                return targets.contains(LayoutOptionData.Target.PORTS);
            } else if (e instanceof ElkLabel) {
                return targets.contains(LayoutOptionData.Target.LABELS);
            }
        }
        return true;
    };

    /**
     * Whether to clear the layout of each graph element before the new configuration is applied.
     */
    public boolean isClearLayout() {
        return clearLayout;
    }
    
    /**
     * Set whether to clear the layout of each graph element before the new configuration is applied
     * (the default is {@code false}).
     * 
     * @return {@code this}
     */
    public LayoutConfigurator setClearLayout(final boolean doClearLayout) {
        this.clearLayout = doClearLayout;
        return this;
    }
    
    /**
     * Adds a filter that is queried for each combination of graph elements and options. An option
     * value is applied only if the filter matches. If no filter is set, all values are applied.
     * 
     * @return {@code this}
     */
    public LayoutConfigurator addFilter(final IOptionFilter filter) {
        this.optionFilters.add(filter);
        return this;
    }
    
    /**
     * Returns the list of filters that have been added via {@link #addFilter(IOptionFilter)} or have been inherited
     * from another {@link LayoutConfigurator} via {@link #overrideWith(LayoutConfigurator)}.
     */
    protected List<IOptionFilter> getFilters() {
        return optionFilters;
    }
    
    /**
     * Add and return a property holder for the given element. If such a property holder is
     * already present, the previous instance is returned.
     */
    public IPropertyHolder configure(final ElkGraphElement element) {
        MapPropertyHolder result = elementOptionMap.get(element);
        if (result == null) {
            result = new MapPropertyHolder();
            elementOptionMap.put(element, result);
        }
        return result;
    }
    
    /**
     * Return the stored property holder for the given element, or {@code null} if none is present.
     */
    public IPropertyHolder getProperties(final ElkGraphElement element) {
        return elementOptionMap.get(element);
    }
    
    /**
     * Add and return a property holder for the given element class. If such a property holder is
     * already present, the previous instance is returned.
     * 
     * It is possible to configure options for every graph element by using {@code KGraphElement.class}
     * as argument. Such options are overridden if a more specific type is configured as well,
     * for instance {@code KNode.class}.
     */
    public IPropertyHolder configure(final Class<? extends ElkGraphElement> elementClass) {
        MapPropertyHolder result = classOptionMap.get(elementClass);
        if (result == null) {
            result = new MapPropertyHolder();
            classOptionMap.put(elementClass, result);
        }
        return result;
    }
    
    /**
     * Return the stored property holder for the given element class, or {@code null} if none is present.
     */
    public IPropertyHolder getProperties(final Class<? extends ElkGraphElement> elementClass) {
        return classOptionMap.get(elementClass);
    }

    /**
     * Apply this layout configurator to the given graph element.
     */
    @Override
    public void visit(final ElkGraphElement element) {
        if (clearLayout) {
            element.getProperties().clear();
        }
        MapPropertyHolder combined = findClassOptions(element);
        // implicitly overwrite options specified for a class with options specified for the specific element
        combined.copyProperties(getProperties(element));
        applyProperties(element, combined);
    }
    
    /**
     * Apply all properties held in {@code properties} to {@code element}.
     */
    @SuppressWarnings("unchecked")
    protected void applyProperties(final ElkGraphElement element, final IPropertyHolder properties) {
        if (properties != null) {
            List<IOptionFilter> filters = getFilters();
            for (Map.Entry<IProperty<?>, Object> entry : properties.getAllProperties().entrySet()) {
                boolean accept = filters.stream()
                        .allMatch(filter -> filter.accept(element, entry.getKey()));
                if (accept) {
                    Object value = entry.getValue();
                    if (value instanceof Cloneable) {
                        Object clone = ElkReflect.clone(value);
                        if (clone != null) {
                            value = clone;
                        }
                    }
                    element.setProperty((IProperty<Object>) entry.getKey(), value);
                }
            }
        }
    }
    
    /**
     * To allow the configuration of layout options using interfaces and super-interfaces, we have to manually check the
     * type hierarchy when visiting a graph element. Since KGraph's type hierarchy is quite small, we do this case by
     * case. The order of the cases is important, since configurations for more specific types should override the
     * general case.
     * 
     * @return the most specific {@link MapPropertyHolder} fitting the passed {@code element}'s type.
     */
    private MapPropertyHolder findClassOptions(final ElkGraphElement element) {
        MapPropertyHolder combined = new MapPropertyHolder();
        
        // order is important as it reflects inheritance 
        
        if (element instanceof ElkGraphElement) {
            combined.copyProperties(classOptionMap.get(ElkGraphElement.class));
        }
        
        if (element instanceof ElkShape) {
            combined.copyProperties(classOptionMap.get(ElkShape.class));
        }
        
        if (element instanceof ElkLabel) {
            combined.copyProperties(classOptionMap.get(ElkLabel.class));
            // cannot be anything of the following, we can return
            return combined;
        }
        
        if (element instanceof ElkConnectableShape) {
            combined.copyProperties(classOptionMap.get(ElkConnectableShape.class));
        }
        
        if (element instanceof ElkNode) {
            combined.copyProperties(classOptionMap.get(ElkNode.class));
            // cannot be anything of the following, we can return
            return combined;
        }
        
        if (element instanceof ElkPort) {
            combined.copyProperties(classOptionMap.get(ElkPort.class));
            // cannot be anything of the following, we can return
            return combined;
        }
        
        if (element instanceof ElkEdge) {
            combined.copyProperties(classOptionMap.get(ElkEdge.class));
        }
        
        return combined;
    }
    
    /**
     * Copy all options from the given configurator to this one, possibly overriding the own options. The
     * {@link IOptionFilter}s of this configurator are cleared and filled with the filters of {@code other}.
     * 
     * @return {@code this}
     */
    public LayoutConfigurator overrideWith(final LayoutConfigurator other) {
        for (Map.Entry<ElkGraphElement, MapPropertyHolder> entry : other.elementOptionMap.entrySet()) {
            MapPropertyHolder thisHolder = this.elementOptionMap.get(entry.getKey());
            if (thisHolder == null) {
                thisHolder = new MapPropertyHolder();
                this.elementOptionMap.put(entry.getKey(), thisHolder);
            }
            thisHolder.copyProperties(entry.getValue());
        }
        for (Map.Entry<Class<? extends ElkGraphElement>, MapPropertyHolder> entry : other.classOptionMap.entrySet()) {
            MapPropertyHolder thisHolder = this.classOptionMap.get(entry.getKey());
            if (thisHolder == null) {
                thisHolder = new MapPropertyHolder();
                this.classOptionMap.put(entry.getKey(), thisHolder);
            }
            thisHolder.copyProperties(entry.getValue());
        }
        this.clearLayout = other.clearLayout;
        this.optionFilters.clear();
        this.optionFilters.addAll(other.optionFilters);
        return this;
    }

}
