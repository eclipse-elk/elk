/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui.views;


import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.ui.LayoutOptionLabelProvider;
import org.eclipse.elk.core.ui.Messages;
import org.eclipse.elk.core.validation.GraphIssue;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.jface.viewers.CellEditor;
import org.eclipse.jface.viewers.ComboBoxCellEditor;
import org.eclipse.jface.viewers.ICellEditorValidator;
import org.eclipse.jface.viewers.ILabelProvider;
import org.eclipse.jface.viewers.TextCellEditor;
import org.eclipse.swt.SWT;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.ui.views.properties.IPropertyDescriptor;
import org.eclipse.ui.views.properties.IPropertySheetEntry;

/**
 * A property descriptor for layout options.
 * 
 * <p>Here's a small peculiarity concerning the descriptions of layout options: calling
 * {@link #getDescription()} only returns the first sentence of an option's description,
 * to be displayed in the status bar. To retrieve the full description, call
 * {@link #getLongDescription()}.</p>
 *
 * @author msp
 */
public class LayoutPropertyDescriptor implements IPropertyDescriptor {

    /** The property source from which this descriptor was created. */
    private final LayoutPropertySource source;
    /** the layout option data associated with this property descriptor. */
    private final LayoutOptionData optionData;
    /** option targets applicable to the currently selected diagram element. */
    private final Set<LayoutOptionData.Target> elementTargets;
    /** the label provider for this property descriptor. */
    private LayoutOptionLabelProvider labelProvider;
    
    /** There's another option that allows to select a layout algorithm other than {@link CoreOptions#ALGORITHM}. 
     *  To avoid a dependency to that plugin, the option is hard-coded here. */
    public static final String DISCO_LAYOUT_ALG_ID =
            "org.eclipse.elk.disco.componentCompaction.componentLayoutAlgorithm";
    /**
     * Creates a layout property descriptor based on a specific layout option.
     * 
     * @param source the property source from which this descriptor was created
     * @param theoptionData the layout option data
     * @param theelementTargets option targets applicable to the currently selected diagram element
     */
    public LayoutPropertyDescriptor(final LayoutPropertySource source,
            final LayoutOptionData theoptionData,
            final Set<LayoutOptionData.Target> theelementTargets) {
        this.source = source;
        this.optionData = theoptionData;
        this.elementTargets = theelementTargets;
    }
    
    /**
     * @return the layout option data associated with this property descriptor
     */
    public LayoutOptionData getOptionData() {
        return optionData;
    }
    
    @Override
    public CellEditor createPropertyEditor(final Composite parent) {
        switch (optionData.getType()) {
        case STRING:
            if (CoreOptions.ALGORITHM.equals(optionData)) {
                return new AlgorithmCellEditor(parent, CoreOptions.ALGORITHM.getId());
            } else if (optionData.getId().equals(DISCO_LAYOUT_ALG_ID)) {
                return new AlgorithmCellEditor(parent, DISCO_LAYOUT_ALG_ID);
            } else {
                return new TextCellEditor(parent);
            }
        case INT:
            CellEditor intEditor = new TextCellEditor(parent);
            intEditor.setValidator(new ICellEditorValidator() {
                public String isValid(final Object value) {
                    try {
                        int i = Integer.parseInt((String) value);
                        List<GraphIssue> issues = source.validatePropertyValue(optionData, i);
                        if (issues != null && !issues.isEmpty()) {
                            return issues.get(0).getMessage();
                        }
                        return null;
                    } catch (NumberFormatException exception) {
                        return Messages.getString("elk.ui.6");
                    }
                }
            });
            return intEditor;
        case DOUBLE:
            CellEditor floatEditor = new TextCellEditor(parent);
            floatEditor.setValidator(new ICellEditorValidator() {
                public String isValid(final Object value) {
                    try {
                        double d = Double.parseDouble((String) value);
                        List<GraphIssue> issues = source.validatePropertyValue(optionData, d);
                        if (issues != null && !issues.isEmpty()) {
                            return issues.get(0).getMessage();
                        }
                        return null;
                    } catch (NumberFormatException exception) {
                        return Messages.getString("elk.ui.7");
                    }
                }
            });
            return floatEditor;
        case BOOLEAN:
            return new ComboBoxCellEditor(parent, optionData.getChoices(), SWT.READ_ONLY);
        case ENUM:
            return new ComboBoxCellEditor(parent, getEnumerationChoices(optionData), SWT.READ_ONLY);
        case ENUMSET:
            return new MultipleOptionsCellEditor(parent, getEnumerationChoices(optionData), true);
        case OBJECT:
            return new TextCellEditor(parent);
        default:
            return null;
        }
    }
    
    /**
     * Returns a string array representing the enumeration choices of the given layout option. The choices are
     * decorated with suffixes if they are advanced, experimental, or both.
     * 
     * @param data the layout option to generate choices for.
     * @return array of decorated choices.
     */
    private String[] getEnumerationChoices(final LayoutOptionData data) {
        String[] choices = new String[data.getEnumValueCount()];
        
        for (int i = 0; i < choices.length; i++) {
            Enum<?> value = data.getEnumValue(i);
            
            // Check if the value is advanced or experimental
            boolean isAdvanced = ElkGraphUtil.isAdvancedPropertyValue(value);
            boolean isExperimental = ElkGraphUtil.isExperimentalPropertyValue(value);
            
            if (isAdvanced && isExperimental) {
                choices[i] = value.toString() + " (advanced and experimental)";
            } else if (isAdvanced) {
                choices[i] = value.toString() + " (advanced)";
            } else if (isExperimental) {
                choices[i] = value.toString() + " (experimental)";
            } else {
                choices[i] = value.toString();
            }
        }
        
        return choices;
    }

    @Override
    public String getCategory() {
        // Compute the intersection of the targets registered for the layout option
        // and those determined for the graph element currently in focus.
        Set<LayoutOptionData.Target> targets = EnumSet.copyOf(optionData.getTargets());
        if (elementTargets != null) {
            targets.retainAll(elementTargets);
        }
        // In normal cases, the intersection should contain exactly one element.
        if (targets.size() == 1) {
            switch (targets.iterator().next()) {
            case PARENTS:
                if (elementTargets != null && elementTargets.contains(LayoutOptionData.Target.NODES)) {
                    return Messages.getString("elk.ui.74");
                }
                return Messages.getString("elk.ui.73");
            case NODES:
                return Messages.getString("elk.ui.75");
            case EDGES:
                return Messages.getString("elk.ui.76");
            case PORTS:
                return Messages.getString("elk.ui.77");
            case LABELS:
                return Messages.getString("elk.ui.78");
            }
        }
        return null;
    }

    @Override
    public String getDescription() {
        // We only return the first sentence of the description. The longer description
        // can be retrieved by calling getLongDescription.
        String description = optionData.getDescription();
        int firstPeriodIndex = description.indexOf('.');
        
        if (firstPeriodIndex > 0) {
            return description.substring(0, firstPeriodIndex + 1);
        } else {
            return description;
        }
    }
    
    /**
     * Returns the option's full description. This is in contrast to {@link #getDescription()},
     * which only returns the description's first sentence.
     * 
     * @return the full description.
     */
    public String getLongDescription() {
        return optionData.getDescription();
    }

    @Override
    public String getDisplayName() {
        return optionData.getName();
    }

    @Override
    public String[] getFilterFlags() {
        if (optionData.getVisibility() == LayoutOptionData.Visibility.ADVANCED) {
            return new String[] { IPropertySheetEntry.FILTER_ID_EXPERT };
        }
        return null;
    }

    @Override
    public Object getHelpContextIds() {
        return null;
    }

    @Override
    public Object getId() {
        return optionData.getId();
    }

    @Override
    public ILabelProvider getLabelProvider() {
        if (labelProvider == null) {
            labelProvider = new LayoutOptionLabelProvider(optionData);
        }
        return labelProvider;
    }

    @Override
    public boolean isCompatibleWith(final IPropertyDescriptor anotherProperty) {
        return anotherProperty instanceof LayoutPropertyDescriptor
                && this.optionData.getId().equals(anotherProperty.getId());
    }

}
