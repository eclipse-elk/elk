/*******************************************************************************
 * Copyright (c) 2009, 2017 Kiel University and others.
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
package org.eclipse.elk.core.data;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.core.util.IDataObject;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;
import org.eclipse.elk.graph.util.ElkReflect;

import com.google.common.collect.Lists;

/**
 * Data type used to store information for a layout option.
 * 
 * @author msp
 */
public final class LayoutOptionData implements ILayoutMetaData, IProperty<Object>,
        Comparable<IProperty<?>> {

    /** Enumeration of data types for layout options. */
    public enum Type {
        /** undefined type. */
        UNDEFINED,
        /** boolean type. */
        BOOLEAN,
        /** integer type. */
        INT,
        /** string type. */
        STRING,
        /** double type. */
        DOUBLE,
        /** enumeration type. */
        ENUM,
        /** enumeration set type. */
        ENUMSET,
        /** {@link IDataObject} type. */
        OBJECT;
    }

    /** Enumeration of target elements for layout options. */
    public enum Target {
        /** parents target (hierarchical nodes). */
        PARENTS,
        /** nodes target. */
        NODES,
        /** edges target. */
        EDGES,
        /** ports target. */
        PORTS,
        /** labels target. */
        LABELS;
    }
    
    /** Enumeration of visibility options for layout options. */
    public enum Visibility {
        /** The option shall always be visible in the UI. */
        VISIBLE,
        /** The option shall be visible only for advanced users. */
        ADVANCED,
        /** The option shall never be visible in the UI. */
        HIDDEN;
    }

    /** identifier of the layout option. */
    private final String id;
    /** the group this layout option is associated with. Note that the group is included in the {@code id}. */
    private final String group;
    /** legacy identifiers of this option. */
    private final String[] legacyIds;
    /** the default value of this option. */
    private final Object defaultValue;
    /** the class that represents this option type. */
    private final Class<?> clazz;
    /** type of the layout option. */
    private final Type type;
    /** user friendly name of the layout option. */
    private final String name;
    /** a description to be displayed in the UI. */
    private final String description;
    /** configured targets. */
    private final Set<Target> targets;
    /** dependencies to other layout options. */
    private final List<Pair<LayoutOptionData, Object>> dependencies = Lists.newLinkedList();
    /** cached value of the available choices. */
    private String[] choices;
    /** visibility in the UI. */
    private final Visibility visibility;
    /** the lower bound for option values. */
    private Object lowerBound;
    /** the upper bound for option values. */
    private Object upperBound;
    
    /**
     * Create a layout option data entry.
     */
    private LayoutOptionData(final Builder builder) {
        this.id = builder.id;
        this.group = builder.group;
        this.name = builder.name;
        this.description = builder.description;
        this.defaultValue = builder.defaultValue;
        this.lowerBound = builder.lowerBound;
        this.upperBound = builder.upperBound;
        this.type = builder.type;
        this.clazz = builder.clazz;
        if (builder.targets == null) {
            this.targets = EnumSet.noneOf(Target.class);
        } else {
            this.targets = builder.targets;
        }
        this.visibility = builder.visibility;
        this.legacyIds = builder.legacyIds;
    }
    
    /**
     * Checks whether the enumeration class is set correctly. This method must
     * not be called for options other than of type {@code enum} and {@code enumset}.
     */
    private void checkEnumClass() {
        if (clazz == null || !clazz.isEnum()) {
            throw new IllegalStateException("Enumeration class expected for layout option " + id);
        }
    }

    /**
     * Checks whether the {@link IDataType} class is set correctly and creates an instance.
     * This method must not be called for options other than of type 'object'.
     * 
     * @return an instance of the data object
     */
    private IDataObject createDataInstance() {
        if (clazz == null 
                // elkjs-exclude-start
                // If this is a problem, it should be observed in the Java world, omit the check for GWT
                || !IDataObject.class.isAssignableFrom(clazz)
                // elkjs-exclude-end
                ) {
            throw new IllegalStateException("IDataType class expected for layout option " + id);
        }
        Object instance = ElkReflect.newInstance(clazz);
        if (instance == null) {
            throw new IllegalStateException("Couldn't create new instance of property '" + id + "'. "
                    + "Make sure its type is registered with the " + ElkReflect.class.getSimpleName()
                    + " utility class.");
        }
        return (IDataObject) instance;
    }
    
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof LayoutOptionData) {
            return this.id.equals(((LayoutOptionData) obj).id);
        } else if (obj instanceof IProperty<?>) {
            return this.id.equals(((IProperty<?>) obj).getId());
        } else {
            return false;
        }
    }
    
    @Override
    public int hashCode() {
        return id.hashCode();
    }

    @Override
    public int compareTo(final IProperty<?> other) {
        return id.compareTo(other.getId());
    }

     @Override
    public String toString() {
        return "Layout Option: " + id;
    }
     
     /**
      * Checks whether or not the {@link #parseValue(String)} will be conceivably able to parse Strings into
      * on option value conforming to this option data's data type. This is the case if {@link #getType()}
      * returns something other than {@link Type#OBJECT} and {@link Type#UNDEFINED} or if {@link #getType()}
      * returns {@link Type#OBJECT} and {@link #getClass()} returns something which implements {@link IDataObject}.
      * 
      * @return {@code true} if String-represented values of this option can be parsed.
      */
     public boolean canParseValue() {
         if (type == Type.UNDEFINED) {
             return false;
         } else if (type == Type.OBJECT) {
             return clazz != null 
                     // elkjs-exclude-start
                     // If this is a problem, it should be observed in the Java world, omit the check for GWT  
                     && IDataObject.class.isAssignableFrom(clazz)
                     // elkjs-exclude-end
                     ; // SUPPRESS CHECKSTYLE Whitespace
         } else {
             return true;
         }
     }

    /**
     * Parses a string value for this layout option.
     * 
     * @param valueString a serialized value
     * @return an instance of the corresponding correctly typed value, or
     *         {@code null} if the given value string is invalid
     */
    @SuppressWarnings({ "unchecked", "rawtypes" })
    public Object parseValue(final String valueString) {
        if (valueString == null || valueString.equals("null")) {
            return null;
        }
        
        // if this open data instance is an ENUMSET, we need to allow empty strings
        // to denote that no enumeration value is selected. Otherwise, we forbid empty strings
        if (valueString.length() == 0 && type != Type.ENUMSET) {
            return null;
        }
        
        switch (type) {
        case BOOLEAN:
            if (valueString.equalsIgnoreCase("true")) {
                return Boolean.TRUE;
            } else if (valueString.equalsIgnoreCase("false")) {
                return Boolean.FALSE;
            } else {
                return null;
            }
        case INT:
            try {
                return Integer.valueOf(valueString);
            } catch (NumberFormatException exception) {
                return null;
            }
        case DOUBLE:
            try {
                return Double.valueOf(valueString);
            } catch (NumberFormatException exception) {
                return null;
            }
        case STRING:
            return valueString;
        case ENUM:
            checkEnumClass();
            return enumForString(valueString);
        case ENUMSET:
            checkEnumClass();
            return enumSetForStringArray((Class<? extends Enum>) clazz, valueString);
        case OBJECT:
            try {
                IDataObject value = createDataInstance();
                value.parse(valueString);
                return value;
            } catch (IllegalArgumentException exception) {
                return null;
            }
        default:
            throw new IllegalStateException("Invalid type set for this layout option.");
        }
    }
    
    /**
     * Tries to turn the given string representation into a set over the enumeration of
     * the type given by the class type parameter. The parameter is supposed to be the
     * {@link #clazz} attribute, but has to be given here explicitly for type inference
     * reasons. The string consists of multiple parts, with each part following the
     * convention specified in the comment of {@link #enumForString(String)}. The format
     * of the string is something like {@code [a, b, c]}.
     * 
     * @param leClazz the enumeration class.
     * @param leString the string to convert.
     * @return the enumeration set.
     */
    @SuppressWarnings({ "rawtypes", "unchecked" })
    private <E extends Enum<E>> EnumSet<E> enumSetForStringArray(final Class<E> leClazz,
            final String leString) {
        
        // raw type for gwt compilation 
        //  (avoid a reflections cast in the 'set.add(o)' below)
        EnumSet set = EnumSet.noneOf(leClazz);
        
        // break the value string into its different components and iterate over them;
        // the string will be of the form "[a, b, c]"
        String[] components = leString.split("[\\[\\]\\s,]+");
        for (String component : components) {
            // Check for empty strings
            if (component.trim().length() == 0) {
                continue;
            }
            
            Object o = enumForString(component);
            
            if (o == null) {
                // we were unable to get the enumeration instance, so the whole enumset
                // is pointless; let's all take a moment of silence to mourn
                return null;
            } else {
                // add the enumeration object to the set
                set.add(o);
            }
        }
        
        return set;
    }

    /**
     * Tries to turn the given string representation into an enumeration of the type given by
     * the {@link #clazz} attribute. The string may either specify one of the named items of
     * the enumeration, or a number referring to one of the enumeration's items by its ordinal
     * number.
     * 
     * @param leString the string to convert to an enumeration.
     * @return the enumeration.
     */
    @SuppressWarnings("unchecked")
    private Object enumForString(final String leString) {
        try {
            @SuppressWarnings("rawtypes")
            Enum<?> value = Enum.valueOf((Class<? extends Enum>) clazz, leString);
            return value;
        } catch (IllegalArgumentException exception) {
            // the value could not be parsed as enumeration constant, try as integer
            try {
                int index = Integer.parseInt(leString);
                Object[] constants = clazz.getEnumConstants();
                if (index >= 0 && index < constants.length) {
                    return constants[index];
                }
            } catch (NumberFormatException e) {
                // ignore exception and return null
            }
            return null;
        }
    }
    
    /**
     * Creates a default-default value for this layout option which is used if no default value was
     * specified for a layout option. In contrast to {@link #getDefault()}, this never returns
     * {@code null} for options with type other than 'object'.
     * 
     * @return a default-default value, depending on the option type
     */
    @SuppressWarnings({ "unchecked", "rawtypes" })
    public Object getDefaultDefault() {
        switch (type) {
        case STRING:
            return "";
        case BOOLEAN:
            return Boolean.FALSE;
        case INT:
            return Integer.valueOf(0);
        case DOUBLE:
            return Double.valueOf(0.0);
        case ENUM:
            checkEnumClass();
            Enum<?>[] enums = ((Class<Enum>) clazz).getEnumConstants();
            return enums[0];
        case ENUMSET:
            checkEnumClass();
            return EnumSet.noneOf(((Class<Enum>) clazz));
        case OBJECT:
            return null;
        default:
            throw new IllegalStateException("Invalid type set for this layout option.");
        }
    }
    
    /** choices for boolean type options. */
    private static final String[] BOOLEAN_CHOICES = { "false", "true" };

    /**
     * Creates an array of choices that can be selected by the user to set a
     * value for this option. This makes only sense for enumeration type
     * or boolean type options.
     * 
     * @return an array of values to be displayed for the user
     */
    @SuppressWarnings({ "unchecked", "rawtypes" })
    public String[] getChoices() {
        if (choices == null) {
            switch (type) {
            case ENUM:
            case ENUMSET:
                checkEnumClass();
                Enum<?>[] enums = ((Class<Enum>) clazz).getEnumConstants();
                choices = new String[enums.length];
                for (int i = 0; i < enums.length; i++) {
                    choices[i] = enums[i].toString();
                }
                break;
            case BOOLEAN:
                choices = BOOLEAN_CHOICES;
                break;
            default:
                choices = new String[0];
            }
        }
        return choices;
    }
    
    /**
     * Returns the number of items in the enumeration class. Calling this method is only meaningful if the option type
     * is either {@link Type#ENUM} or {@link Type#ENUMSET}.
     * 
     * @return the number of enumeration constants or 0 if the option type is not enumeration related.
     */
    public int getEnumValueCount() {
        switch (type) {
        case ENUM:
        case ENUMSET:
            checkEnumClass();
            @SuppressWarnings({ "unchecked", "rawtypes" })
            int count = ((Class<? extends Enum>) clazz).getEnumConstants().length;
            return count;
        default:
            return 0;
        }
    }

    /**
     * Returns the enumeration value for a given index. Calling this method is only meaningful if the option type
     * is either {@link Type#ENUM} or {@link Type#ENUMSET}.
     * 
     * @param intValue zero-based index of the enumeration value
     * @return the corresponding enumeration value, or {@code null} if the option type is not enumeration related.
     */
    public Enum<?> getEnumValue(final int intValue) {
        switch (type) {
        case ENUM:
        case ENUMSET:
            checkEnumClass();
            @SuppressWarnings({ "unchecked", "rawtypes" })
            Enum<?>[] enums = ((Class<? extends Enum>) clazz).getEnumConstants();
            return enums[intValue]; 
        default:
            return null;
        }
    }

    /**
     * Returns the set of layout option targets, which determine the types of graph elements
     * the option can be related to.
     * 
     * @return the layout option targets
     */
    public Set<Target> getTargets() {
        return targets;
    }
    
    /**
     * Returns the dependencies to other layout options. The option should only be made visible if
     * at least one of the dependencies is met.
     * 
     * <p>The returned list should not be modified except in {@link ILayoutMetaDataProvider} implementations.</p>
     * 
     * @return the options from which this option depends together with the expected values
     */
    public List<Pair<LayoutOptionData, Object>> getDependencies() {
        return dependencies;
    }

    @Override
    public String getId() {
        return id;
    }
    
    /**
     * @return the group
     */
    public String getGroup() {
        return group;
    }

    /**
     * Returns the type.
     *
     * @return the type
     */
    public Type getType() {
        return type;
    }

    /**
     * Returns the name.
     *
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * Returns the description.
     *
     * @return the description
     */
    public String getDescription() {
        return description;
    }
    
    /**
     * Returns the default value of this layout option.
     * 
     * @return the default value.
     */
    public Object getDefault() {
        // Clone the default value if it's a Cloneable. One would have to use reflection for this to work
        // properly (classes implementing Cloneable are not required to make their clone() method
        // public, so we need to check if they have such a method and invoke it via reflection, which
        // results in ugly and unchecked type casting)
        // HOWEVER, since GWT doesn't support reflection, a helper class is used that is able to clone 
        // those objects that are viable as property values
        if (defaultValue instanceof Cloneable) {
            Object clone = ElkReflect.clone(defaultValue);
            if (clone == null) {
                throw new IllegalStateException(
                        "Couldn't clone property '" + id + "'. " + "Make sure it's type is registered with the "
                                + ElkReflect.class.getSimpleName() + " utility class.");
            }
            return clone;
        } else {
            return defaultValue;
        }
    }

    @Override
    @SuppressWarnings("unchecked")
    public Comparable<? super Object> getLowerBound() {
        if (lowerBound instanceof Comparable<?>) {
            return (Comparable<Object>) lowerBound;
        }
        return (Comparable<Object>) Property.NEGATIVE_INFINITY;
    }

    /**
     * Sets the lower bound for layout option values.
     * 
     * @param lowerBound the lowerBound to set
     */
    public void setLowerBound(final Object lowerBound) {
        this.lowerBound = lowerBound;
    }

    @Override
    @SuppressWarnings("unchecked")
    public Comparable<? super Object> getUpperBound() {
        if (upperBound instanceof Comparable<?>) {
            return (Comparable<Object>) upperBound;
        }
        return (Comparable<Object>) Property.POSITIVE_INFINITY;
    }

    /**
     * Sets the upper bound for layout option values.
     * 
     * @param upperBound the upperBound to set
     */
    public void setUpperBound(final Object upperBound) {
        this.upperBound = upperBound;
    }
    
    /**
     * Returns the option type class.
     * 
     * @return the type class
     */
    public Class<?> getOptionClass() {
        return clazz;
    }

    /**
     * Returns the visibility of this option in the UI.
     * 
     * @return the visibility setting
     */
    public Visibility getVisibility() {
        return visibility;
    }

    /**
     * @return the legacyIds, may be {@code null}
     */
    public String[] getLegacyIds() {
        return legacyIds;
    }
    
    /**
     * Builder for {@link LayoutOptionData} instances.
     */
    public static class Builder {
        
        private String id;
        private String group;
        private String[] legacyIds;
        private Object defaultValue;
        private Class<?> clazz;
        private Type type;
        private String name;
        private String description;
        private Set<Target> targets;
        private Visibility visibility;
        private Object lowerBound;
        private Object upperBound;
        
        /**
         * Create an instance with the configured values.
         */
        public LayoutOptionData create() {
            return new LayoutOptionData(this);
        }
        
        /**
         * Configure the {@link LayoutOptionData#getId() id}.
         */
        public Builder id(final String aid) {
            this.id = aid;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getGroup() group}.
         */
        public Builder group(final String agroup) {
            this.group = agroup;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getLegacyIds() legacyIds}.
         */
        public Builder legacyIds(final String... alegacyIds) {
            this.legacyIds = alegacyIds;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getDefault() defaultValue}.
         */
        public Builder defaultValue(final Object adefaultValue) {
            this.defaultValue = adefaultValue;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getOptionClass() optionClass}.
         */
        public Builder optionClass(final Class<?> aclazz) {
            this.clazz = aclazz;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getType() type}.
         */
        public Builder type(final Type atype) {
            this.type = atype;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getName() name}.
         */
        public Builder name(final String aname) {
            this.name = aname;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getDescription() description}.
         */
        public Builder description(final String adescription) {
            this.description = adescription;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getTargets() targets}.
         */
        public Builder targets(final Set<Target> atargets) {
            this.targets = atargets;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getVisibility() visibility}.
         */
        public Builder visibility(final Visibility avisibility) {
            this.visibility = avisibility;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getLowerBound() lowerBound}.
         */
        public Builder lowerBound(final Object alowerBound) {
            this.lowerBound = alowerBound;
            return this;
        }
        
        /**
         * Configure the {@link LayoutOptionData#getUpperBound() upperBound}.
         */
        public Builder upperBound(final Object aupperBound) {
            this.upperBound = aupperBound;
            return this;
        }
        
    }
    
}
