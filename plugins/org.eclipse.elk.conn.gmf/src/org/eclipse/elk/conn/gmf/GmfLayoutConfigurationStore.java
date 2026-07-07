/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.conn.gmf;

import java.util.Collection;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.data.LayoutOptionData.Target;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.service.ILayoutConfigurationStore;
import org.eclipse.emf.edit.domain.EditingDomain;
import org.eclipse.gef.EditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.AbstractBorderItemEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.CompartmentEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.ConnectionEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.DiagramEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.DiagramRootEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.IGraphicalEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.LabelEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.NoteEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.ShapeNodeEditPart;
import org.eclipse.gmf.runtime.notation.NotationFactory;
import org.eclipse.gmf.runtime.notation.StringValueStyle;
import org.eclipse.gmf.runtime.notation.View;
import org.eclipse.ui.IWorkbenchPart;

import com.google.common.collect.Maps;
import com.google.inject.Inject;

/**
 * A layout configuration that stores layout options in the notation model of GMF diagrams.
 *
 * @author msp
 */
public class GmfLayoutConfigurationStore implements ILayoutConfigurationStore {
    
    /**
     * Provider for GMF layout configuration stores.
     */
    public static final class Provider implements ILayoutConfigurationStore.Provider {
        
        @Inject
        private IEditPartFilter editPartFilter;

        @Override
        public ILayoutConfigurationStore get(final IWorkbenchPart workbenchPart, final Object context) {
            if (context instanceof EditPart) {
                try {
                    return new GmfLayoutConfigurationStore((EditPart) context, editPartFilter);
                } catch (IllegalArgumentException e) {
                    // Fall back to null
                }
            }
            return null;
        }
        
    }
    
    /** Prefix for all layout options. */
    public static final String PREFIX = "layout:";
    
    /** The graphical edit part used as context for this configuration store. */
    private final IGraphicalEditPart editPart;
    /** The edit part filter for excluding specific diagram parts. */
    private final IEditPartFilter editPartFilter;
    
    /**
     * Create a GMF layout configuration store.
     */
    public GmfLayoutConfigurationStore(final EditPart theeditPart, final IEditPartFilter filter) {
        if (theeditPart instanceof CompartmentEditPart) {
            // If the selected object is a compartment, replace it by its parent element
            this.editPart = (IGraphicalEditPart) ((CompartmentEditPart) theeditPart).getParent();
        } else if (theeditPart instanceof IGraphicalEditPart) {
            this.editPart = (IGraphicalEditPart) theeditPart;
        } else if (theeditPart instanceof DiagramRootEditPart) {
            this.editPart = (IGraphicalEditPart) ((DiagramRootEditPart) theeditPart).getContents();
        } else {
            throw new IllegalArgumentException("Not supported: " + theeditPart.toString());
        }
        this.editPartFilter = filter;
    }

    @Override
    public Object getOptionValue(final String optionId) {
        View view = editPart.getNotationView();
        if (view != null) {
            String key = PREFIX + optionId;
            for (Object obj : view.getStyles()) {
                if (obj instanceof StringValueStyle) {
                    StringValueStyle style = (StringValueStyle) obj;
                    if (key.equals(style.getName())) {
                        String result = style.getStringValue();
                        if (result != null) {
                            return result;
                        }
                    }
                }
            }
        }
        Map<String, Object> defaultOptions = getDefaultOptions();
        if (defaultOptions != null) {
            return defaultOptions.get(optionId);
        }
        return null;
    }
    
    /**
     * Return a map of default options. The base implementation only assigns the option
     * {@link CoreOptions#COMMENT_BOX} to instances of {@link NoteEditPart}. More defaults can
     * be specified in subclasses.
     */
    protected Map<String, Object> getDefaultOptions() {
        Map<String, Object> result = null;
        if (editPart instanceof NoteEditPart) {
            result = Maps.newHashMapWithExpectedSize(1);
            result.put(CoreOptions.COMMENT_BOX.getId(), true);
        }
        return result;
    }

    @Override
    @SuppressWarnings("unchecked")
    public void setOptionValue(final String optionId, final String value) {
        View view = editPart.getNotationView();
        if (view != null) {
            String key = PREFIX + optionId;
            if (value == null) {
                removeValue(key, view);
            } else {
                for (Object obj : view.getStyles()) {
                    if (obj instanceof StringValueStyle) {
                        StringValueStyle style = (StringValueStyle) obj;
                        if (key.equals(style.getName())) {
                            style.setStringValue(value);
                            return;
                        }
                    }
                }
        
                StringValueStyle style = NotationFactory.eINSTANCE.createStringValueStyle();
                style.setName(key);
                style.setStringValue(value);
                view.getStyles().add(style);
            }
        }
    }

    /**
     * Remove an option from the given notation view.
     */
    protected void removeValue(final String key, final View view) {
        Iterator<?> iter = view.getStyles().iterator();
        while (iter.hasNext()) {
            Object obj = iter.next();
            if (obj instanceof StringValueStyle
                    && key.equals(((StringValueStyle) obj).getName())) {
                iter.remove();
            }
        }
    }
    
    @Override
    public Collection<String> getAffectedOptions() {
        Set<String> options = new HashSet<String>();
        View view = editPart.getNotationView();
        if (view != null) {
            for (Object obj : view.getStyles()) {
                if (obj instanceof StringValueStyle) {
                    StringValueStyle style = (StringValueStyle) obj;
                    String key = style.getName();
                    if (key != null && key.startsWith(PREFIX)) {
                        options.add(key.substring(PREFIX.length()));
                    }
                }
            }
        }
        Map<String, Object> defaultOptions = getDefaultOptions();
        if (defaultOptions != null) {
            options.addAll(defaultOptions.keySet());
        }
        return options;
    }
    
    @Override
    public EditingDomain getEditingDomain() {
        return editPart.getEditingDomain();
    }

    @Override
    public Set<Target> getOptionTargets() {
        if (editPart instanceof AbstractBorderItemEditPart) {
            // This is a border item, i.e. a port 
            return EnumSet.of(LayoutOptionData.Target.PORTS);
        } else if (editPart instanceof ShapeNodeEditPart) {
            // This is a node
            Set<LayoutOptionData.Target> partTarget = EnumSet.of(LayoutOptionData.Target.NODES);
            // Check whether the node is a parent
            if (findContainingEditPart(editPart) != null) {
                partTarget.add(LayoutOptionData.Target.PARENTS);
            }
            return partTarget;
        } else if (editPart instanceof ConnectionEditPart) {
            // This is a connection, i.e. an edge
            return EnumSet.of(LayoutOptionData.Target.EDGES);
        } else if (editPart instanceof LabelEditPart) {
            // This is a label
            return EnumSet.of(LayoutOptionData.Target.LABELS);
        } else if (editPart instanceof DiagramEditPart) {
            // This is a diagram
            return EnumSet.of(LayoutOptionData.Target.PARENTS);
        }
        return EnumSet.noneOf(LayoutOptionData.Target.class);
    }
    
    /**
     * Finds the edit part that contains layoutable children, if there are any. The returned
     * edit part is either the parent edit part itself or one of its compartments.
     * 
     * @param nodeEditPart a node edit part
     * @return the edit part that contains other node edit parts, or {@code null} if there is none
     */
    protected IGraphicalEditPart findContainingEditPart(final IGraphicalEditPart nodeEditPart) {
        for (Object child : nodeEditPart.getChildren()) {
            if (editPartFilter.filter((EditPart) child)) {
                if (child instanceof ShapeNodeEditPart) {
                    return nodeEditPart;
                } else if (child instanceof CompartmentEditPart) {
                    for (Object grandChild : ((CompartmentEditPart) child).getChildren()) {
                        if (grandChild instanceof ShapeNodeEditPart
                                && editPartFilter.filter((EditPart) grandChild)) {
                            return (IGraphicalEditPart) child;
                        }
                    }
                }
            }
        }
        return null;
    }
    
    @Override
    public ILayoutConfigurationStore getParent() {
        EditPart container = getContainer();
        if (container != null) {
            return new GmfLayoutConfigurationStore(container, editPartFilter);
        }
        return null;
    }
    
    /**
     * Determines the container edit part representing a parent node of the context edit part.
     */
    protected EditPart getContainer() {
        if (editPart instanceof AbstractBorderItemEditPart) {
            // The container is the parent of the port's containing node
            return (IGraphicalEditPart) editPart.getParent().getParent();
        } else if (editPart instanceof ShapeNodeEditPart) {
            // This is a node
            return (IGraphicalEditPart) editPart.getParent();
        } else if (editPart instanceof ConnectionEditPart) {
            EditPart sourcePart = ((ConnectionEditPart) editPart).getSource();
            if (sourcePart instanceof AbstractBorderItemEditPart) {
                // The source element is a port
                return sourcePart.getParent().getParent();
            } else if (sourcePart != null) {
                // The source element is a node
                return sourcePart.getParent();
            }
        } else if (editPart instanceof LabelEditPart) {
            // This is a label
            IGraphicalEditPart containerEditPart = (IGraphicalEditPart) editPart.getParent();
            if (containerEditPart instanceof ConnectionEditPart) {
                // We have an edge label, so apply the same container rule as for edges
                EditPart sourcePart = ((ConnectionEditPart) containerEditPart).getSource();
                if (sourcePart instanceof AbstractBorderItemEditPart) {
                    return (IGraphicalEditPart) sourcePart.getParent().getParent();
                } else if (sourcePart != null) {
                    return (IGraphicalEditPart) sourcePart.getParent();
                }
            } else if (containerEditPart instanceof AbstractBorderItemEditPart) {
                // We have a port label, so apply the same container rule as for ports
                return (IGraphicalEditPart) containerEditPart.getParent().getParent();
            } else if (containerEditPart instanceof ShapeNodeEditPart) {
                // We have a node label
                return (IGraphicalEditPart) containerEditPart.getParent();
            }
            return containerEditPart;
        }
        return null;
    }
    
}
